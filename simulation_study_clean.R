library(DPpack)
library(glmnet)
###################################################################
### Functions for simulating data and training/evaluating model ###
###################################################################

#' Generate simulation data
#'
#' @param N Number of simulated examples.
#' @param p Number of predictors.
#' @param propensities Vector of probabilities of individual being randomly
#'   assigned to -1 or 1, respectively.
#'   
#' @return List consisting of
#'   X: N x p feature data frame consisting of independently drawn features 
#'     from unif(0, 1)
#'   A: N-length vector of randomized treatment assignments
#'   B: N-length vector of treatment benefits
#'   f: Ground truth decision function
#'   upper.bounds: p-length vector of upper bounds on feature columns
#'   lower.bounds: p-length vector of lower bounds on feature columns
#'   
#' @examples 
#' N <- 300 # Number of training datapoints
#' p <- 10 # Number of features
#' training.data <- generate.data(N, p)
generate.data <- function(N, p, propensities=c(0.5, 0.5)){
  # Generate X
  X <- matrix(NaN, nrow=N, ncol=p)
  for (i in 1:N){
    X[i, ] <- runif(p)
  }
  X <- data.frame(X)
  
  ## NOTE: Upper and lower bounds on features should ideally be selected 
  ##  independent of the data. Ideally, there would be some natural lower and 
  ##  upper limit to each feature. For our simulations, the features are 
  ##  simulated from the unif(0, 1) distribution, so the limits are as stated 
  ##  below.
  upper.bounds <- rep(1, p)
  lower.bounds <- rep(0, p)
  
  # Generate A
  A <- sample(c(-1, 1), size=N, replace=TRUE, prob=propensities)
  
  # Ground truth function
  f <- function(x){
    if (is.null(dim(x))){
      return(1 + x[1] + x[2] - 1.8*x[3] - 2.2*x[4])
    } else{
      return(1 + x[,1] + x[,2] - 1.8*x[,3] - 2.2*x[,4])
    }
  }
  
  # Generate B
  mus <- 0.01 + 0.02*X[,4] + 3*A*f(X)
  std.dev <- 0.5
  B <- numeric(N)
  for (i in 1:N){
    B[i] <- rnorm(1, mean=mus[i], sd=std.dev)
  }
  if (min(B) < 0) B <- B + abs(min(B)) + 0.001
  
  return(list(X=X, A=A, B=B, f=f, upper.bounds=upper.bounds, 
              lower.bounds=lower.bounds))
}

#' Train wSVM model using simulated training data
#'
#' @param training.data List of training data information in format of output of
#'   generate.data function.
#' @param eps Privacy budget epsilon.
#' @param gamma Regularization constant.
#' @param kernel Kernel for SVM. Can be 'linear' for linear kernel or 'Gaussian'
#'   for Gaussian (radial) kernel.
#' @param D Dimension of Gaussian kernel approximation. Only used if kernel is
#'   'Gaussian'.
#' @param weights.upper.bound Upper bound for weights necessary for DP.
#' @param propensity Single number corresponding to larger of probabilities of
#'   assigning individual to A=1 or A=-1.
#' @param add.bias Boolean indicating whether to include a bias term in the
#'   resulting predictor function.
#'   
#' @return Fitted OWL model.
#' 
#' @examples 
#'   N <- 300
#'   p <- 10
#'   training.data <- generate.data(N, p)
#'   eps <- 1
#'   gamma <- 50
#'   kernel <- 'linear'
#'   D <- 15 # unused if kernel is not 'Gaussian'
#'   wub <- 30
#'   owl.model <- train.model(training.data, eps, gamma, kernel, D, wub)
train.model <- function(training.data, eps, gamma, kernel, D, 
                        weights.upper.bound, propensity=0.5, add.bias=TRUE){
  # Redefine target variable to align with svmDP requirements that y = 0 or 1
  y <- training.data$A
  y[y<0] <- 0
  # Bounds on training data (0 and 1 since drawn from U(0, 1))
  p <- ncol(training.data$X)
  # Weights from B and propensity score
  weights <- training.data$B/propensity
  
  # Construct object
  wsvm <- svmDP$new("l2", eps, gamma, perturbation.method="output", kernel, 
                    D=D)
  
  # Train model
  wsvm$fit(training.data$X, y, training.data$upper.bounds, 
           training.data$lower.bounds, add.bias, weights, weights.upper.bound)
  
  return(wsvm)
}

#' Helper function for estimating empirical value function on validation data
#'    (best to call evaluate.model directly, which calls this function)
#'
#' @param trained.model Trained OWL model of form of output of train.model
#'   function.
#' @param val.X Feature data frame from validation data.
#' @param val.A Treatment assignment vector from validation data.
#' @param val.B Treatment benefit vector from validation data.
#' @param propensity Single number corresponding to larger of probabilities of
#'   assigning individual to A=1 or A=-1.
#' @param add.bias Boolean indicating whether to include a bias term in the
#'   predictor function.
#'   
#' @return Estimated empirical value function.
#'   
#' @examples 
#'   trained.model <- train.model(...) # This is output of train.model function
#'   validation.data <- generate.data(500, 10) # Validation data
#'   estimate.value.function(trained.model, validation.data$X, 
#'                           validation.data$A, validation.data$B, 0.5, TRUE)
estimate.value.function <- function(trained.model, val.X, val.A, val.B,
                                    propensity, add.bias){
  val.A.hat <- trained.model$predict(val.X, add.bias)
  top <- sum((val.A.hat==val.A)*val.B/propensity)/length(val.A.hat)
  bottom <- sum((val.A.hat==val.A)/propensity)/length(val.A.hat)
  return(top/bottom)
}

#' Evaluate fitted OWL model using validation data
#'
#' @param trained.model Trained OWL model of form of output of train.model
#'   function.
#' @param validation.data Validation data of form of output of generate.data
#'   function.
#' @param propensity Single number corresponding to larger of probabilities of
#'   assigning individual to A=1 or A=-1.
#' @param add.bias Boolean indicating whether to include a bias term in the
#'   predictor function.
#' @param ground.truth Boolean indicating whether a ground truth decision 
#'   function is available.
#'
#' @return List containing the percent accuracy of the trained model and the
#'   estimated empirical value function on the validation data.
#'
#' @examples 
#'   trained.model <- train.model(...) # This is output of train.model function
#'   validation.data <- generate.data(500, 10) # Validation data
#'   results <- evaluate.model(trained.model, validation.data)
#'   results$percent.accuracy # Gets percent accuracy of model on validation set
#'   results$val.func # Gets estimated value function on validation set
evaluate.model <- function(trained.model, validation.data, propensity=0.5,
                           add.bias=TRUE, ground.truth=TRUE){
  if (ground.truth){
    # Compare assigned treatment level to ground truth best treatment label
    validation.A.hat <- trained.model$predict(validation.data$X, add.bias=add.bias)
    # Convert to -1 and 1 for comparison
    validation.A.hat[validation.A.hat==0] <- -1 
    
    # Ground truth best treatment option
    validation.A.best <- sign(validation.data$f(validation.data$X))
    
    percent.accuracy <- 100*sum(validation.A.hat==validation.A.best)/
      length(validation.A.hat)
  } else percent.accuracy <- "N/A"
  
  # Empirical value function
  val.func <- estimate.value.function(trained.model, validation.data$X,
                                      validation.data$A, validation.data$B,
                                      propensity, add.bias)
  
  return(list(percent.accuracy=percent.accuracy, val.func=val.func))
}

##############################################
### Feature selection via LASSO and/or PCA ###
##############################################

## NOTE: The number of features seems to greatly impact the performance of the 
##  model, as more features requires a larger vector of noise and therefore more
##  noise overall. Thus, we ideally want to use as few features as possible. We
##  suggest doing this through a combination of LASSO regression and PCA. The
##  suggested workflow would be
##  1) If there are few features to begin with (~5 or so), try training the 
##    model using all features.
##  2) If there are more features, attempt to use weighted LASSO regression 
##    (example below) to select only a few important features to use. Using the
##    original features can maintain the explanability benefits.
##  3) If it appears that many features are significant, use PCA to first 
##    transform features to significant components, then use weighted LASSO 
##    regression to select only a few important ones. In the simulated datasets,
##    we assumed the features were independent, so we don't use PCA here.
##  The feature selection will not be released, so it can be done without DP.

# Generate sample data for feature selection with 10 features
propensities <- c(0.5, 0.5)
fs.data <- generate.data(500, 10, propensities)

# # Direct tuning with access to ground truth function (only possible for 
# #   simulated data)
# cv.mod <- cv.glmnet(as.matrix(fs.data$X),
#                     fs.data$f(fs.data$X)) # family='binomial' also works well
# coef(cv.mod, s = "lambda.1se")

# Indirect tuning with only access to observed benefits and random assignments
cv.mod <- cv.glmnet(as.matrix(fs.data$X), fs.data$A,
                    weights=fs.data$B/propensities[1], family='binomial')
coef(cv.mod, s="lambda.1se")

# Based on results, we should only use first 4 predictors, which matches the
#   ground truth function

################################
### Hyperparameter selection ###
################################

## NOTE: The privacy budget eps should be chosen.
eps <- 5

## NOTE: The upper bound on the weights should ideally be selected independent
##  of the data. Ideally, there would be some natural lower and upper limit to
##  the benefit metric. For DP, the benefit values should be first shifted so
##  that they are all nonnegative, then an upper bound for these shifted values
##  should be set. If the upper bound is set too tight (i.e. lower than the 
##  true upper bound), the benefits that are too large will simply be clipped to
##  the upper bound. Selecting this value as tight possible to the true largest
##  benefit value in the dataset produces best results. The upper bound for the
##  weights is then the benefit upper bound divided by the larger of the 
##  probabilities of being assigned to -1 or 1 for the treatment. For the
##  P(A=-1) = P(A=1) = 0.5 case, the weight upper bound should be 2*benefit 
##  upper bound.
wub <- 30 # Global/public weights upper bound (sup{benefit/propensity})

## NOTE: The regularization constant has an effect on the performance of the
##  OWL algorithm and needs to be carefully chosen to optimize performance.
##  With a large enough dataset, a portion of the privacy budget can be devoted
##  to selecting this value. However, with a dataset of the size we are likely 
##  to see in the OWL setting, this would likely be difficult. The best 
##  regularization constant depends on the value of the privacy budget eps
##  and the sample size. It may also depend on the dataset, but that is not yet
##  clear. As a first attempt, I recommend using the 
##  `Regularization constant data.csv` table to select a value for gamma based 
##  on the closest possible values of eps and N to the real data.
gamma <- 60 # Value associated with eps=5 and N=500

## NOTE: It is also necessary to select lower/upper bounds for feature vectors 
##  x. For the simulated data example here, data are simulated with natural 
##  lower/upper bounds of 0 and 1. See generate.data() function for more info.
##  Real-world data will likely not satisfy this and thus different bound need
##  to be set. If the lower/upper bounds are set too tight (e.g., lower than 
##  the true upper bound), the bounds that are too large/small will simply be 
##  clipped to the given bound. Selecting this value as tight possible to the 
##  true smallest/largest feature value in the dataset produces best results.
# lower.bounds <- 
# upper.bounds <- 

## NOTE: All of our experiments were with a 'linear' kernel, which worked well 
##  in the simulations. It may also be worthwhile to try the 'Gaussian' kernel,
##  in which case you will also need to select the approximation dimensionality 
##  D.
kernel <- 'linear'
# kernel <- 'Gaussian'
D <- 15 # Unused if kernel=='linear'

################################
### Example run with dataset ###
################################

#### NOTE: Replace this part with your data after feature selection. To easily
####  apply the functions in this script, be sure your data is in the same
####  format as this example with the exception of the ground truth function.
####  The format should be a list with 
####    $X N x p data.frame of features
####    $A N-dimensional vector of treatment assignments
####    $B N-dimensional vector of treatment benefits, shifted so B_i>=0
####    $upper.bounds p-dimensional vector of upper bounds on features
####    $lower.bounds p-dimensional vector of lower bounds on features
N.train <- 300
p <- 4 # Chosen from feature selection
training.data <- generate.data(N.train, p)

N.val <- 200
validation.data <- generate.data(N.val, p)
####
####

owl.model <- train.model(training.data = training.data,
                         eps = eps,
                         gamma = gamma,
                         kernel = kernel,
                         D = D,
                         weights.upper.bound = wub)

results <- evaluate.model(trained.model = owl.model, 
                          validation.data = validation.data,
                          # NOTE: Set to FALSE if ground truth unknown
                          ground.truth = TRUE) 

print(paste("Epsilon:", eps))
print(paste("Prediction Accuracy:", results$percent.accuracy, "%"))
print(paste("Estimated Value Function:", results$val.func))

