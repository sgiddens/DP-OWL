library(ggplot2)
library(ggbreak)
library(collections)
library(DPpack)
library(glmnet)
###########################
### Simulation Settings ###
###########################
set.seed(100) # Set this for consistency in generating figures for DP-OWL paper

benefit.scale <- 1
wub <- 30/benefit.scale # Global/public upper bound on weights (sup{benefit/propensity})
normal.X <- FALSE
trim.B <- FALSE
use.01 <- FALSE
trim.B.bounds <- c(4, 10)
if (trim.B) wub <- (trim.B.bounds[2]-trim.B.bounds[1])/0.5 # assuming propensity=0.5
if (trim.B) use.01 <- FALSE
# wub <- 6 # For when X is normalized
# normal.X <- TRUE
# set.seed(42)

# Partition data and use parallel composition to aggregate results by 
#   averaging coefficients
n.divisions <- 1

# Subset data and use privacy amplification to aggregate results by 
#   averaging coefficients 
# NOTE: n.divisions above also indicates the number of subsets
subset.rate <- 1

# Use SVM without weighting 
unweighted.svm <- FALSE

# Use negative weights
neg.weights <- FALSE

if (neg.weights) wub <- 15

#####################
### Generate data ###
#####################
generate.data <- function(N, p, propensities=c(0.5, 0.5), normalize.X=FALSE){
  # Generate X
  X <- matrix(NaN, nrow=N, ncol=p)
  for (i in 1:N){
    tmp <- runif(p)
    if (normalize.X) tmp <- tmp/sqrt(p+1) # Normalize X
    X[i, ] <- tmp
  }
  X <- data.frame(X)

  # Generate A
  A <- sample(c(-1, 1), size=N, replace=TRUE, prob=propensities)
  if (use.01) A <- sample(c(0, 1), size=N, replace=TRUE, prob=propensities) 

  # Ground truth function
  f <- function(x){
    if (is.null(dim(x))){
      const <- 1
      if (normalize.X) const <- const/sqrt(p+1)
      out <- const + x[1] + x[2] - 1.8*x[3] - 2.2*x[4]
      return(out)
    } else{
      const <- 1
      if (normalize.X) const <- const/sqrt(p+1)
      out <- const + x[,1] + x[,2] - 1.8*x[,3] - 2.2*x[,4]
      return(out)
    }
  }

  # Generate B
  mus <- 0.01 + 0.02*X[,4] + 3*A*f(X)
  std.dev <- 0.5
  B <- numeric(N)
  for (i in 1:N){
    B[i] <- rnorm(1, mean=mus[i], sd=std.dev)
  }
  if (!neg.weights){
    if (min(B) < 0) B <- B + abs(min(B)) + 0.001
  }
  
  if (trim.B){
    B.mask.small <- (B<trim.B.bounds[1])
    B[B.mask.small] <- trim.B.bounds[1]
    B.mask.large <- (B>trim.B.bounds[2])
    B[B.mask.large] <- trim.B.bounds[2]
    B <- B - trim.B.bounds[1]
  }
  
  # # Discretize B
  # B[B<2.5] <- 0.001
  # B[2.5<=B & B<5] <- 0.25
  # B[5<=B & B<10] <- 0.5
  # B[10<=B & B<12.5] <- 0.75
  # B[12.5<=B] <- 1
  
  # Normalize B to reduce sensitivity - relative benefit matters (not absolute)
  B <- B/benefit.scale

  return(list(X=X, A=A, B=B, f=f))
}

### Single example ###
# N <- 300 # Number of training datapoints
# p <- 10 # Number of features
# training.data <- generate.data(N, p)

###################################
### Feature Selection via LASSO ###
###################################
# ## NOTE: Ideally, one would use PCA on a real dataset prior to feature selection
# 
# # NOTE: This feature selection will not be released, so can be done without DP
# # Generate sample data for feature selection
# propensities <- c(0.5, 0.5)
# fs.data <- generate.data(1000, 10, propensities)
# 
# # Direct tuning with access to ground truth function
# cv.mod <- cv.glmnet(as.matrix(fs.data$X),
#                     fs.data$f(fs.data$X)) # family='binomial' also works well
# coef(cv.mod, s = "lambda.1se")
# 
# cv.mod <- cv.glmnet(as.matrix(fs.data$X),
#                     fs.data$B)
# coef(cv.mod, s = "lambda.1se")
# 
# # Indirect tuning with only access to observed benefits and random assignments
# 
# cv.mod <- cv.glmnet(as.matrix(fs.data$X), fs.data$A,
#                     weights=fs.data$B/propensities[1], family='binomial')
# coef(cv.mod, s="lambda.1se")
# 
# # Indicates to use only p=4

################################
### Train wSVM model with DP ###
################################
train.model <- function(training.data, eps, gamma, kernel, D,
                        weights.upper.bound, propensity=0.5, add.bias=TRUE,
                        n.divisions=1, subset.rate=1.0){
  # Redefine target variable to align with svmDP requirements
  y <- training.data$A
  y[y<0] <- 0
  # Bounds on training data (0 and 1 since drawn from U(0, 1))
  p <- ncol(training.data$X)
  upper.bounds <- rep(1, p)
  if (normal.X) upper.bounds <- upper.bounds/sqrt(p+1)
  lower.bounds <- rep(0, p)
  # Weights from B and propensity score
  weights <- training.data$B/propensity
  
  # Shuffle all indices
  all.indices <- sample(1:length(training.data$A))
  
  # Partition data
  division.indices <- split(all.indices, 1:n.divisions)
  
  # Subset data
  if (subset.rate<1){
    division.indices <- list()
    for (i in 1:n.divisions){
      division.indices[[i]] <- all.indices[
        runif(length(training.data$A))<=subset.rate]
    }
    
    eps <- eps/(n.divisions*subset.rate)
  }
  
  # All coefficients
  all.coeffs <- matrix(nrow=n.divisions, ncol=ncol(training.data$X)+1)
  if (!add.bias) all.coeffs <- matrix(nrow=n.divisions, 
                                      ncol=ncol(training.data$X))
  for (i in 1:n.divisions){
    idxs <- division.indices[[i]]
    
    # set.seed(NULL)
    # Construct object
    wsvm <- svmDP$new("l2", eps, gamma, perturbation.method="output", kernel, 
                      D=D)
    
    # Train model
    if (unweighted.svm){
      wsvm$fit(training.data$X[idxs,], y[idxs], upper.bounds, lower.bounds, 
               add.bias)
    } else{
      wsvm$fit(training.data$X[idxs,], y[idxs], upper.bounds, lower.bounds, 
               add.bias, weights[idxs], weights.upper.bound)
    }
    all.coeffs[i,] <- wsvm$coeff
  }
  out.wsvm.model <- svmDP$new("l2", eps, gamma, perturbation.method="output", 
                              kernel, D=D)
  if (unweighted.svm){
    out.wsvm.model$fit(training.data$X, y, upper.bounds, lower.bounds,
                       add.bias)
  } else{
    out.wsvm.model$fit(training.data$X, y, upper.bounds, lower.bounds,
                       add.bias, weights, weights.upper.bound)
  }
  out.wsvm.model$coeff <- apply(all.coeffs, 2, mean)
  return(out.wsvm.model)
}

### Single example ###
# eps <- 1 # NOTE: Can use without DP by setting to Inf
# gamma <- 10
# # kernel <- "Gaussian" # Radial kernal
# kernel <- "linear"
# D <- 15 # Dimension of projected kernel approximation space
# # This needs to be a previously fixed public value
# #  Theoretical bounds from the simulation settings show the weights are upper
# #  bounded by approximately 35. This number can be reduced if clipping weights
# #  on the higher end is acceptable
# weights.upper.bound <- wub # This essentially should be sup{training.data$B/propensity}
# wsvm <- train.model(training.data, eps, gamma, kernel, D, weights.upper.bound)

##############################
### Benefit clipping stuff ###
##############################

# ### Histogram to determine benefit cutoffs ###
# ### NOTE: This method of doing it is not correct
# # a <- generate.data(1000, 10)
# # correct.A <- as.numeric(a$f(a$X)>0)
# # d <- data.frame(B=a$B, true.A=correct.A)
# # ggplot(data=d, aes(x=B, color=as.factor(true.A), fill=as.factor(true.A))) +
# #   geom_histogram(position='dodge')
# 
# # Platt scaling to get probability of assigning A=1 vs B plot
# train.data <- generate.data(1000, 10)
# test.data <- generate.data(5000, 10)
# # Treatment assignment model
# svm.model <- train.model(train.data, Inf, 1, "linear", 5, wub)
# 
# y <- factor(sign(train.data$f(train.data$X)))
# x <- svm.model$predict(train.data$X, TRUE, TRUE)
# platt.scaling <- glm(y ~ x, family = "binomial")
# 
# # plot(train.data$B, predict(platt.scaling, train.data$X, type="response"))
# plot.df <- data.frame(Benefit=train.data$B, 
#                       PA1=predict(platt.scaling, train.data$X, type="response"))
# ggplot(plot.df,
#        aes(Benefit, PA1)) + geom_point() + ylab("P(A=1)") +
#   theme(text = element_text(size = 16))
# # Results of this plot show a good clipping range is 4<B<10

###########################
### Evaluate wSVM model ###
###########################
# Estimate value function
estimate.value.function <- function(trained.model, val.X, val.A, val.B,
                                    propensity, add.bias){
  val.A.hat <- trained.model$predict(val.X, add.bias)
  if (!use.01) val.A.hat[val.A.hat==0] <- -1
  top <- sum((val.A.hat==val.A)*val.B/propensity)/length(val.A.hat)
  bottom <- sum((val.A.hat==val.A)/propensity)/length(val.A.hat)
  return(top/bottom)
}

evaluate.model <- function(trained.model, validation.data, propensity=0.5,
                           add.bias=TRUE){
  # Compare assigned treatment level to ground truth best treatment label
  validation.A.hat <- trained.model$predict(validation.data$X, add.bias=add.bias)
  # Convert to -1 and 1 for comparison
  validation.A.hat[validation.A.hat==0] <- -1 

  # Ground truth best treatment option
  validation.A.best <- sign(validation.data$f(validation.data$X))

  percent.accuracy <- 100*sum(validation.A.hat==validation.A.best)/
    length(validation.A.hat)

  # Empirical value function
  val.func <- estimate.value.function(trained.model, validation.data$X,
                                      validation.data$A, validation.data$B,
                                      propensity, add.bias)

  return(list(percent.accuracy=percent.accuracy, val.func=val.func))
}

### Single example ###
# validation.data <- generate.data(5000, p)
# results <- evaluate.model(wsvm, validation.data)
#
# print(paste("Epsilon:", eps))
# print(paste("Prediction Accuracy:", results$percent.accuracy, "%"))
# print(paste("Estimated Value Function:", results$val.func))

###################################################################
### Tune Regularization Constants for each Training Sample Size ###
###################################################################
# NOTE: The regularization constant gamma will be tuned using an independent 
#   public dataset with the same training size. In practice, with access to some 
#   public dataset that is believed to come from a similar distribution as the 
#   private one, we would either undersample (if public dataset bigger) or 
#   bootstrap (if public dataset smaller) the public dataset to get a dataset of 
#   the same size as the publicly known size of the private training dataset.
#   We would then tune the value using that sample size and use the tuned value 
#   to train the private dataset.

# This function returns a dataframe with each regularization constant in the 
#   tuning set and their computed validation accuracies and estimated empirical
#   value functions for selection. Best to choose regularization constant
#   based on estimated empirical value function as that does not assume 
#   knowledge of true treatment assignment function f.
evaluate.reg.constant <- function(tuning.set, n.simulations, public.data, 
                                  sample.size, eps, kernel="linear", D=15, 
                                  weights.upper.bound=wub, propensity=0.5, 
                                  add.bias=TRUE, train.proportion=0.8){
  # Determine with/without replacement (over/undersampling)
  public.data.size <- length(public.data$A)
  if (floor(train.proportion*public.data.size) < sample.size){
    replace <- TRUE
  } else replace <- FALSE
  
  mean.accs <- numeric(length(tuning.set))
  sd.accs <- numeric(length(tuning.set))
  mean.val.funcs <- numeric(length(tuning.set))
  sd.val.funcs <- numeric(length(tuning.set))
  for (i in 1:length(tuning.set)){
    gamma <- tuning.set[i]
    print(paste("Evaluating gamma =", gamma, "..."))
    tmp.accs <- numeric(n.simulations)
    tmp.val.funcs <- numeric(n.simulations)
    for (j in 1:n.simulations){
      # Get training/validation indices
      train.indices <- sample(floor(train.proportion*public.data.size), 
                              size=sample.size,
                              replace=replace)
      validation.indices <- sample((floor(
        train.proportion*public.data.size)+1):public.data.size, 
        size=floor((1-train.proportion)*public.data.size))
      
      # Get training/validation data
      training.data <- list()
      training.data$X <- public.data$X[train.indices, ]
      training.data$A <- public.data$A[train.indices]
      training.data$B <- public.data$B[train.indices]
      
      validation.data <- list()
      validation.data$X <- public.data$X[validation.indices, ]
      validation.data$A <- public.data$A[validation.indices]
      validation.data$B <- public.data$B[validation.indices]
      validation.data$f <- public.data$f
      
      trained.model <- train.model(training.data, eps, gamma, kernel, D,
                                   weights.upper.bound, propensity, add.bias)
      results <- evaluate.model(trained.model, validation.data, propensity,
                                add.bias)
      tmp.accs[j] <- results$percent.accuracy
      tmp.val.funcs[j] <- results$val.func
    }
    mean.accs[i] <- mean(tmp.accs)
    sd.accs[i] <- sd(tmp.accs)
    mean.val.funcs[i] <- mean(tmp.val.funcs)
    sd.val.funcs[i] <- sd(tmp.val.funcs)
  }
  data.frame(Regularization.Constant=tuning.set,
             Mean.Percent.Validation.Accuracy=mean.accs,
             Std.Dev.Percent.Validation.Accuracy=sd.accs,
             Mean.Estimated.Validation.Value=mean.val.funcs,
             Std.Dev.Estimated.Validation.Value=sd.val.funcs)
}
public.data <- generate.data(1000, 10)
# public.data <- fs.data
# public.data$X <- public.data$X[,1:4]

get.reg.constant <- function(tuned.reg.constants, eps, sample.size){
  eps.mask <- tuned.reg.constants$eps==eps
  tmp1 <- tuned.reg.constants[eps.mask,]
  sample.size.mask <- tmp1$sample.size==sample.size
  tmp2 <- tmp1[sample.size.mask,]
  out <- tmp2$gamma
  if (length(out)!=1) stop(paste("data.frame of tuned regularization constants",
              " does not yield adequate output for these inputs"))
  out
}

#####
###### Tuning for eps = 0.1
{
### Sample size 200 (eps=0.1, p=10) ###
sample.size <- 200
eps <- 0.1
tuned.reg.constants <- data.frame(eps=eps, sample.size=sample.size,
                          gamma=50) # 50 determined by below parameter tuning
#####
# gammas <- 10*(1:10) - 9 # Not much utility of note anywhere at this eps
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=2/10)

#####
### Sample size 500 (eps=0.1, p=10) ###
sample.size <- 500
eps <- 0.1
if (neg.weights) gam <- 125 # p=4
else gam <- 75
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                          sample.size=sample.size,
                          gamma=gam)) # gam determined by below parameter tuning
#####
# # gammas <- 10*(1:10) - 9
# gammas <- 10*(3:18)
# # gammas <- 10*(1:11)+40 # Not much utility of note anywhere at this eps
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=0.5)

#####
### Sample size 800 (eps=0.1, p=10) ###
sample.size <- 800
eps <- 0.1
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                          sample.size=sample.size,
                          gamma=75)) # 75 determined by below parameter tuning
#####
# # gammas <- 5*(1:20)+15 # Not much utility anywhere at this eps
# gammas <- 5*(1:20) - 4
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=0.5)

#####
### Sample size 1000 (eps=0.1, p=10) ###
sample.size <- 1000
eps <- 0.1
if (neg.weights) gam <- 175 # p=4
else gam <- 110
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                          sample.size=sample.size,
                          gamma=gam)) # gam determined by below parameter tuning
#####
# # gammas <- 5*(1:20)+35 # Not much utility anywhere at this eps
# # gammas <- 5*(1:40) - 4
# # gammas <- 10*(1:25)
# gammas <- 10*(17:30)
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=0.5)
#####
### Sample size 1500 (eps=0.1, p=10) ###
sample.size <- 1500
eps <- 0.1
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                          sample.size=sample.size,
                          gamma=110)) # 110 determined by below parameter tuning
#####
# # gammas <- 5*(1:20)+35 # Not much utility anywhere at this eps
# gammas <- 5*(1:25) - 4
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=0.5)
#####
### Sample size 2000 (eps=0.1, p=10) ###
sample.size <- 2000
eps <- 0.1
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                           sample.size=sample.size,
                           gamma=175)) # 175 determined by below parameter tuning
#####
# gammas <- 5*(1:40)
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=0.5)
#####
### Sample size 2500 (eps=0.1, p=10) ###
sample.size <- 2500
eps <- 0.1
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                            sample.size=sample.size,
                            gamma=225)) # 225 determined by below parameter tuning
#####
# gammas <- 5*(1:40)+70
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=0.5)
}

#####
###### Tuning for eps = 0.4
{
### Sample size 50 (eps=0.4, p=10) ###
sample.size <- 50
eps <- 0.4
gam <- 90
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                         sample.size=sample.size,
                         gamma=gam)) # determined by below parameter tuning
#####
# gammas <- 5*(1:50)
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=2/10)
}

#####
###### Tuning for eps = 0.5
{
### Sample size 200 (eps=0.5, p=10) ###
sample.size <- 200
eps <- 0.5
gam <- 90
if (trim.B) gam <- 90
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                            sample.size=sample.size,
                            gamma=gam)) # determined by below parameter tuning
#####
# # gammas <- 5*(1:20)+35
# # gammas <- 5*(1:20)+75
# gammas <- 5*(1:25)
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=2/10)

#####
### Sample size 500 (eps=0.5, p=10) ###
sample.size <- 500
eps <- 0.5
if (neg.weights) gam <- 100 # p=4
else gam <- 150
if (trim.B) gam <- 50
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                          sample.size=sample.size,
                          gamma=gam)) # determined by below parameter tuning
#####
# # gammas <- 5*(1:20)+35
# # gammas <- 5*(1:20)+75
# gammas <- 5*(1:25)
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=0.5)

#####
### Sample size 800 (eps=0.5, p=10) ###
sample.size <- 800
eps <- 0.5
gam <- 145
if (trim.B) gam <- 50
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                          sample.size=sample.size,
                          gamma=gam)) # determined by below parameter tuning
#####
# # gammas <- 5*(1:20)+75
# gammas <- 5*(5:35)
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=0.5)

#####
### Sample size 1000 (eps=0.5, p=10) ###
sample.size <- 1000
eps <- 0.5
if (neg.weights) gam <- 150 # p=4
else gam <- 125
if (trim.B) gam <- 100
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                          sample.size=sample.size,
                          gamma=gam)) # determined by below parameter tuning
#####
# # gammas <- 5*(1:20)+75
# gammas <- 5*(1:40)
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=0.5)
#####
### Sample size 1500 (eps=0.5, p=10) ###
sample.size <- 1500
eps <- 0.5
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                        sample.size=sample.size,
                        gamma=125)) # 125 determined by below parameter tuning
#####
# gammas <- 5*(1:20)+75
# # gammas <- 5*(1:25)
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=0.5)
#####
### Sample size 2000 (eps=0.5, p=10) ###
sample.size <- 2000
eps <- 0.5
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                            sample.size=sample.size,
                            gamma=225)) # 225 determined by below parameter tuning
#####
# gammas <- 5*(1:40)+50
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=0.5)
#####
### Sample size 2500 (eps=0.5, p=10) ###
sample.size <- 2500
eps <- 0.5
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                          sample.size=sample.size,
                          gamma=250)) # 250 determined by below parameter tuning
#####
# gammas <- 5*(1:40)+75
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=0.5)
}

#####
###### Tuning for eps = 1
{
### Sample size 200 (eps=1, p=10) ###
sample.size <- 200
eps <- 1
if (neg.weights) gam <- 50 # p=10 and p=4
else gam <- 150
if (trim.B) gam <- 150
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                          sample.size=sample.size,
                          gamma=gam)) # determined by below parameter tuning
#####
# # gammas <- 5*(1:20)+75
# # gammas <- 5*(1:10)+35
# gammas <- 5*(1:40)
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=2/10)

#####
### Sample size 500 (eps=1, p=10) ###
sample.size <- 500
eps <- 1
if (neg.weights){
  gam <- 75
  if (ncol(public.data$X)==4) gam <- 100
} 
else gam <- 115
if (trim.B) gam <- 25
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                          sample.size=sample.size,
                          gamma=gam)) # determined by below parameter tuning
#####
# # gammas <- 5*(1:25)+50
# gammas <- 5*(1:40)
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=0.5)

#####
### Sample size 800 (eps=1, p=10) ###
sample.size <- 800
eps <- 1
gam <- 75
if (ncol(public.data$X)==4) gam <- 125 # p=4
if (trim.B) gam <- 50
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                            sample.size=sample.size,
                            gamma=gam)) # determined by below parameter tuning
#####
# # gammas <- 5*(1:25)+30
# gammas <- 5*(1:40)
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=0.5)

#####
### Sample size 1000 (eps=1, p=10) ###
sample.size <- 1000
eps <- 1
if (neg.weights) {
  gam <- 100
  if (ncol(public.data$X)==4) gam <- 150
}
else gam <- 75
if (trim.B) gam <- 50
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                            sample.size=sample.size,
                            gamma=gam)) # determined by below parameter tuning
#####
# # gammas <- 5*(1:25)+30
# gammas <- 5*(1:40)
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=0.5)
#####
### Sample size 1500 (eps=1, p=10) ###
sample.size <- 1500
eps <- 1
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                            sample.size=sample.size,
                            gamma=125)) # 125 determined by below parameter tuning
#####
# # gammas <- 5*(1:25)+30
# # gammas <- 5*(1:25)
# gammas <- 5*(1:25)+60
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=0.5)
#####
### Sample size 2000 (eps=1, p=10) ###
sample.size <- 2000
eps <- 1
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                          sample.size=sample.size,
                          gamma=225)) # 225 determined by below parameter tuning
#####
# gammas <- 5*(1:40)+100
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=0.5)
#####
### Sample size 2500 (eps=1, p=10) ###
sample.size <- 2500
eps <- 1
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                          sample.size=sample.size,
                          gamma=250)) # 250 determined by below parameter tuning
#####
# gammas <- 5*(1:40)+125
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=0.5)
}

#####
###### Tuning for eps = 2
{
### Sample size 100 (eps=2, p=10) ###
sample.size <- 100
eps <- 2
if (neg.weights) gam <- 40
else gam <- 105
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                            sample.size=sample.size,
                            gamma=gam)) # determined by below parameter tuning
#####
# gammas <- 5*(1:30)
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=2/10)
  
#####
### Sample size 200 (eps=2, p=10) ###
sample.size <- 200
eps <- 2
if (neg.weights) gam <- 30
else gam <- 105
if (trim.B) gam <- 20
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                          sample.size=sample.size,
                          gamma=gam)) # determined by below parameter tuning
#####
# # gammas <- 5*(1:25)+30
# # gammas <- 2*(1:20)
# gammas <- 5*(1:40)
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=2/10)

#####
### Sample size 500 (eps=2, p=10) ###
sample.size <- 500
eps <- 2
gam <- 55
if (neg.weights){
  if (ncol(public.data$X)==4) gam <- 75 # p=4
}
if (trim.B) gam <- 20
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                            sample.size=sample.size,
                            gamma=gam)) # determined by below parameter tuning
#####
# # gammas <- 5*(1:25)+30
# gammas <- 5*(1:25)
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=0.5)

#####
### Sample size 800 (eps=2, p=10) ###
sample.size <- 800
eps <- 2
if (neg.weights) {
  gam <- 50
  if (ncol(public.data$X)==4) gam <- 100
}
else gam <- 100
if (trim.B) gam <- 35
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                          sample.size=sample.size,
                          gamma=gam)) # determined by below parameter tuning
#####
# # gammas <- 5*(1:25)+15
# gammas <- 5*(1:30)
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=0.5)

#####
### Sample size 1000 (eps=2, p=10) ###
sample.size <- 1000
eps <- 2
if (neg.weights) gam <- 100
else gam <- 150
if (trim.B) gam <- 40
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                          sample.size=sample.size,
                          gamma=gam)) # determined by below parameter tuning
#####
# # gammas <- 5*(1:30)+15
# # gammas <- 5*(1:10)+150
# gammas <- 5*(1:30)
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=0.5)
#####
### Sample size 1500 (eps=2, p=10) ###
sample.size <- 1500
eps <- 2
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                            sample.size=sample.size,
                            gamma=125)) # 125 determined by below parameter tuning
#####
# gammas <- 5*(1:40)
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=0.5)
#####
### Sample size 2000 (eps=2, p=10) ###
sample.size <- 2000
eps <- 2
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                          sample.size=sample.size,
                          gamma=225)) # 225 determined by below parameter tuning
#####
# gammas <- 5*(1:40)+100
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=0.5)
#####
### Sample size 2500 (eps=2, p=10) ###
sample.size <- 2500
eps <- 2
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                          sample.size=sample.size,
                          gamma=300)) # 300 determined by below parameter tuning
#####
# gammas <- 5*(1:40)+150
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=0.5)
}

#####
###### Tuning for eps = 4
{
### Sample size 50 (eps=4, p=10) ###
sample.size <- 50
eps <- 4
gam <- 35
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                         sample.size=sample.size,
                         gamma=gam)) # determined by below parameter tuning
#####
# gammas <- 5*(1:50)
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=2/10)
}

#####
###### Tuning for eps = 5
{
### Sample size 200 (eps=5, p=10) ###
sample.size <- 200
eps <- 5
gam <- 30
if (trim.B) gam <- 20
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                            sample.size=sample.size,
                            gamma=gam)) # determined by below parameter tuning
#####
# # gammas <- c(20, 22, 24, 26, 28, 30, 32, 34, 36)
# gammas <- 5*(1:20)
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=2/10)

#####
### Sample size 500 (eps=5, p=10) ###
sample.size <- 500
eps <- 5
if (neg.weights) {
  gam <- 30
  if (ncol(public.data$X)==4) gam <- 60
}
else gam <- 60
if (trim.B) gam <- 20
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                            sample.size=sample.size,
                            gamma=gam)) # determined by below parameter tuning
#####
# # gammas <- c(1, 10, 20, 30, 40, 50, 60, 70) # First attempt - course grid
# # gammas <- c(40, 50, 60, 70, 80, 90, 100, 110) # Next attempt, etc.
# # gammas <- c(50, 55, 60, 65, 70, 75, 80)
# # gammas <- c(50, 52.5, 55, 57.5, 60, 62.5, 65, 67.5, 70)
# gammas <- 5*(1:20)
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=0.5)

#####
### Sample size 800 (eps=5, p=10) ###
sample.size <- 800
eps <- 5
gam <- 50
if (trim.B) gam <- 15
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                            sample.size=sample.size,
                            gamma=gam)) # determined by below parameter tuning
#####
# # gammas <- c(40, 50, 60, 70, 80, 90, 100, 110)
# # gammas <- c(30, 35, 40, 45, 50, 55, 60, 65)
# gammas <- (1:20)*5
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=0.5)

#####
### Sample size 1000 (eps=5, p=10) ###
sample.size <- 1000
eps <- 5
gam <- 75
if (trim.B) gam <- 25
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                            sample.size=sample.size,
                            gamma=gam)) # determined by below parameter tuning
#####
# # gammas <- c(1, 10, 20, 30, 40, 50, 60, 70)
# # gammas <- c(40, 45, 50, 55, 60, 65, 70, 75, 80)
# # gammas <- c(70, 72.5, 75, 77.5, 80, 82.5, 85, 87.5, 90)
# # gammas <- 5*(15:45)
# gammas <- 5*(1:30)
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=0.5)
#####
### Sample size 1500 (eps=5, p=10) ###
sample.size <- 1500
eps <- 5
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                                sample.size=sample.size,
                                gamma=75)) # 75 from below parameter tuning
#####
# gammas <- 5*(1:40)
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=0.5)
#####
### Sample size 2000 (eps=5, p=10) ###
sample.size <- 2000
eps <- 5
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                              sample.size=sample.size,
                              gamma=125)) # 125 from below parameter tuning
#####
# gammas <- 5*(1:50)+50
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=0.5)
#####
### Sample size 2500 (eps=5, p=10) ###
sample.size <- 2500
eps <- 5
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                        sample.size=sample.size,
                        gamma=225)) # 225 from below parameter tuning
#####
# gammas <- 10*(1:35)+50
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=0.5)
}

#####
###### Tuning for eps = 20
{
  ### Sample size 200 (eps=5, p=10) ###
  sample.size <- 200
  eps <- 20
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                        sample.size=sample.size,
                        gamma=20)) # 20 determined by below parameter tuning
  #####
  # gammas <- 5*(1:20)
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=2/10)
  
  #####
  ### Sample size 500 (eps=5, p=10) ###
  sample.size <- 500
  eps <- 20
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                        sample.size=sample.size,
                        gamma=30)) # 30 determined by below parameter tuning
  #####
  # gammas <- 5*(1:20) 
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=0.5)
  
  #####
  ### Sample size 800 (eps=5, p=10) ###
  sample.size <- 800
  eps <- 20
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                  sample.size=sample.size,
                  gamma=50)) # 50 determined by below parameter tuning
  #####
  # gammas <- (1:20)*5
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=0.5)
  
  #####
  ### Sample size 1000 (eps=5, p=10) ###
  sample.size <- 1000
  eps <- 20
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                            sample.size=sample.size,
                            gamma=60)) # 60 unnormalized benefit
  #####
  # gammas <- 5*(1:20)
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=0.5)
  #####
  ### Sample size 1500 (eps=5, p=10) ###
  sample.size <- 1500
  eps <- 20
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                        sample.size=sample.size,
                        gamma=60)) # 60 from below parameter tuning
  #####
  # gammas <- 5*(1:40)
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=0.5)
  #####
  ### Sample size 2000 (eps=5, p=10) ###
  sample.size <- 2000
  eps <- 20
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                          sample.size=sample.size,
                          gamma=100)) # 100 from below parameter tuning
  #####
  # gammas <- 5*(1:40)
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=0.5)
  #####
  ### Sample size 2500 (eps=5, p=10) ###
  sample.size <- 2500
  eps <- 20
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                          sample.size=sample.size,
                          gamma=100)) # 100 from below parameter tuning
  #####
  # gammas <- 10*(1:35)+50
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=0.5)
}

#####
###### Tuning for eps = 50
{
  ### Sample size 200 (eps=5, p=10) ###
  sample.size <- 200
  eps <- 50
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                          sample.size=sample.size,
                          gamma=10)) # 10 determined by below parameter tuning
  #####
  # gammas <- 5*(1:20)
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=2/10)
  
  #####
  ### Sample size 500 (eps=5, p=10) ###
  sample.size <- 500
  eps <- 50
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                          sample.size=sample.size,
                          gamma=30)) # 30 determined by below parameter tuning
  #####
  # gammas <- 5*(1:15)-4
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=0.5)
  
  #####
  ### Sample size 800 (eps=5, p=10) ###
  sample.size <- 800
  eps <- 50
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                        sample.size=sample.size,
                        gamma=30)) # 30 determined by below parameter tuning
  #####
  # gammas <- (1:20)*5
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=0.5)
  
  #####
  ### Sample size 1000 (eps=5, p=10) ###
  sample.size <- 1000
  eps <- 50
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                       sample.size=sample.size,
                       gamma=30)) # 30 determined by below parameter tuning
  #####
  # gammas <- 5*(1:20)
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=0.5)
  #####
  ### Sample size 1500 (eps=5, p=10) ###
  sample.size <- 1500
  eps <- 50
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                     sample.size=sample.size,
                     gamma=30)) # 30 from below parameter tuning
  #####
  # gammas <- 5*(1:40)
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=0.5)
  #####
  ### Sample size 2000 (eps=5, p=10) ###
  sample.size <- 2000
  eps <- 50
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                    sample.size=sample.size,
                    gamma=35)) # 35 from below parameter tuning
  #####
  # gammas <- 5*(1:40)
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=0.5)
  #####
  ### Sample size 2500 (eps=5, p=10) ###
  sample.size <- 2500
  eps <- 50
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                    sample.size=sample.size,
                    gamma=40)) # 40 from below parameter tuning
  #####
  # gammas <- 10*(1:35)+50
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=0.5)
}

#####
###### Tuning for eps = 150
{
  ### Sample size 200 (eps=5, p=10) ###
  sample.size <- 200
  eps <- 150
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                           sample.size=sample.size,
                           gamma=10)) # 10 determined by below parameter tuning
  #####
  # gammas <- (1:20)/2
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=2/10)
  
  #####
  ### Sample size 500 (eps=5, p=10) ###
  sample.size <- 500
  eps <- 150
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                              sample.size=sample.size,
                              gamma=15)) # 15 determined by below parameter tuning
  #####
  # gammas <- 5*(1:15)-4
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=0.5)
  
  #####
  ### Sample size 800 (eps=5, p=10) ###
  sample.size <- 800
  eps <- 150
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                          sample.size=sample.size,
                          gamma=15)) # 15 determined by below parameter tuning
  #####
  # gammas <- (1:20)*5
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=0.5)
  
  #####
  ### Sample size 1000 (eps=5, p=10) ###
  sample.size <- 1000
  eps <- 150
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                          sample.size=sample.size,
                          gamma=15)) # 15 determined by below parameter tuning
  #####
  # gammas <- 5*(1:20)
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=0.5)
  #####
  ### Sample size 1500 (eps=5, p=10) ###
  sample.size <- 1500
  eps <- 150
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                            sample.size=sample.size,
                            gamma=30)) # 30 from below parameter tuning
  #####
  # gammas <- 5*(1:20)
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=0.5)
  #####
  ### Sample size 2000 (eps=5, p=10) ###
  sample.size <- 2000
  eps <- 150
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                        sample.size=sample.size,
                        gamma=20)) # 20 from below parameter tuning
  #####
  # gammas <- 5*(1:20)
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=0.5)
  #####
  ### Sample size 2500 (eps=5, p=10) ###
  sample.size <- 2500
  eps <- 150
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                          sample.size=sample.size,
                          gamma=25)) # 25 from below parameter tuning
  #####
  # gammas <- 5*(1:20)
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=0.5)
}

#####
###### Tuning for eps = 300
{
  ### Sample size 200 (eps=5, p=10) ###
  sample.size <- 200
  eps <- 300
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                          sample.size=sample.size,
                          gamma=10)) # 10 determined by below parameter tuning
  #####
  # gammas <- (1:20)/2
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=2/10)
  
  #####
  ### Sample size 500 (eps=5, p=10) ###
  sample.size <- 500
  eps <- 300
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                          sample.size=sample.size,
                          gamma=10)) # 10 determined by below parameter tuning
  #####
  # gammas <- 5*(1:15)-4
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=0.5)
  
  #####
  ### Sample size 800 (eps=5, p=10) ###
  sample.size <- 800
  eps <- 300
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                          sample.size=sample.size,
                          gamma=10)) # 10 determined by below parameter tuning
  #####
  # gammas <- (1:15)*5-4
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=0.5)
  
  #####
  ### Sample size 1000 (eps=5, p=10) ###
  sample.size <- 1000
  eps <- 300
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                          sample.size=sample.size,
                          gamma=15)) # 15 determined by below parameter tuning
  #####
  # gammas <- 5*(1:15)-4
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=0.5)
  #####
  ### Sample size 1500 (eps=5, p=10) ###
  sample.size <- 1500
  eps <- 300
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                          sample.size=sample.size,
                          gamma=10)) # 10 from below parameter tuning
  #####
  # gammas <- 5*(1:20)
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=0.5)
  #####
  ### Sample size 2000 (eps=5, p=10) ###
  sample.size <- 2000
  eps <- 300
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                       sample.size=sample.size,
                       gamma=40)) # 40 from below parameter tuning
  #####
  # gammas <- 5*(1:15)
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=0.5)
  #####
  ### Sample size 2500 (eps=5, p=10) ###
  sample.size <- 2500
  eps <- 300
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                        sample.size=sample.size,
                        gamma=20)) # 20 from below parameter tuning
  #####
  # gammas <- 5*(1:15)
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=0.5)
}

#####
###### Tuning for eps = 500
{
  ### Sample size 200 (eps=5, p=10) ###
  sample.size <- 200
  eps <- 500
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                  sample.size=sample.size,
                  gamma=5)) # 5 determined by below parameter tuning
  #####
  # gammas <- (1:20)/2
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=2/10)
  
  #####
  ### Sample size 500 (eps=5, p=10) ###
  sample.size <- 500
  eps <- 500
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                   sample.size=sample.size,
                   gamma=5)) # 5 determined by below parameter tuning
  #####
  # gammas <- 5*(1:15)-4
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=0.5)
  
  #####
  ### Sample size 800 (eps=5, p=10) ###
  sample.size <- 800
  eps <- 500
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                          sample.size=sample.size,
                          gamma=5)) # 5 determined by below parameter tuning
  #####
  # gammas <- (1:15)*5-4
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=0.5)
  
  #####
  ### Sample size 1000 (eps=5, p=10) ###
  sample.size <- 1000
  eps <- 500
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                          sample.size=sample.size,
                          gamma=10)) # 10 determined by below parameter tuning
  #####
  # gammas <- 5*(1:15)-4
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=0.5)
  
  #####
  ### Sample size 1500 (eps=5, p=10) ###
  sample.size <- 1500
  eps <- 500
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                            sample.size=sample.size,
                            gamma=15)) # 15 from below parameter tuning
  #####
  # gammas <- 5*(1:20)
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=0.5)
  #####
  ### Sample size 2000 (eps=5, p=10) ###
  sample.size <- 2000
  eps <- 500
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                        sample.size=sample.size,
                        gamma=10)) # 10 from below parameter tuning
  #####
  # gammas <- 5*(1:15)
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=0.5)
  #####
  ### Sample size 2500 (eps=5, p=10) ###
  sample.size <- 2500
  eps <- 500
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                        sample.size=sample.size,
                        gamma=15)) # 15 from below parameter tuning
  #####
  # gammas <- 5*(1:15)
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=0.5)
}

#####
###### Tuning for eps = 800
{
  ### Sample size 200 (eps=5, p=10) ###
  sample.size <- 200
  eps <- 800
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                    sample.size=sample.size,
                    gamma=5)) # 5 determined by below parameter tuning
  #####
  # gammas <- (1:20)/2
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=2/10)
  
  #####
  ### Sample size 500 (eps=5, p=10) ###
  sample.size <- 500
  eps <- 800
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                        sample.size=sample.size,
                        gamma=5)) # 5 determined by below parameter tuning
  #####
  # gammas <- 5*(1:15)-4
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=0.5)
  
  #####
  ### Sample size 800 (eps=5, p=10) ###
  sample.size <- 800
  eps <- 800
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                        sample.size=sample.size,
                        gamma=5)) # 5 determined by below parameter tuning
  #####
  # gammas <- (1:15)*5-4
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=0.5)
  
  #####
  ### Sample size 1000 (eps=5, p=10) ###
  sample.size <- 1000
  eps <- 800
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                        sample.size=sample.size,
                        gamma=5)) # 5 determined by below parameter tuning
  #####
  # gammas <- 5*(1:15)-4
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=0.5)
  
  #####
  ### Sample size 1500 (eps=5, p=10) ###
  sample.size <- 1500
  eps <- 800
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                        sample.size=sample.size,
                        gamma=10)) # 10 from below parameter tuning
  #####
  # gammas <- 5*(1:20)
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=0.5)
  #####
  ### Sample size 2000 (eps=5, p=10) ###
  sample.size <- 2000
  eps <- 800
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                            sample.size=sample.size,
                            gamma=5)) # 5 from below parameter tuning
  #####
  # gammas <- 5*(1:15)-4
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=0.5)
  #####
  ### Sample size 2500 (eps=5, p=10) ###
  sample.size <- 2500
  eps <- 800
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                        sample.size=sample.size,
                        gamma=20)) # 20 from below parameter tuning
  #####
  # gammas <- 5*(1:15)-4
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=0.5)
}

#####
###### Tuning for eps = 1000
{
  ### Sample size 200 (eps=5, p=10) ###
  sample.size <- 200
  eps <- 1000
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                        sample.size=sample.size,
                        gamma=4)) # 4 determined by below parameter tuning
  #####
  # gammas <- 2.5*(1:20)-2
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=2/10)
  
  #####
  ### Sample size 500 (eps=5, p=10) ###
  sample.size <- 500
  eps <- 1000
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                        sample.size=sample.size,
                        gamma=10)) # 10 determined by below parameter tuning
  
  #####
  # # gammas <- 5*(1:20)
  # gammas <- 1:20
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=0.5)
  
  #####
  ### Sample size 800 (eps=5, p=10) ###
  sample.size <- 800
  eps <- 1000
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                          sample.size=sample.size,
                          gamma=1)) # 1 determined by below parameter tuning
  #####
  # gammas <- (1:20)*2.5
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=0.5)
  
  #####
  ### Sample size 1000 (eps=5, p=10) ###
  sample.size <- 1000
  eps <- 1000
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                            sample.size=sample.size,
                            gamma=7)) # 7 determined by below parameter tuning
  #####
  # # gammas <- (1:20)/2
  # gammas <- (1:20)*5
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=0.5)
  #####
  ### Sample size 1500 (eps=5, p=10) ###
  sample.size <- 1500
  eps <- 1000
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                        sample.size=sample.size,
                        gamma=5)) # 5 from below parameter tuning
  #####
  # gammas <- 5*(1:20)-4
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=0.5)
  #####
  ### Sample size 2000 (eps=5, p=10) ###
  sample.size <- 2000
  eps <- 1000
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                          sample.size=sample.size,
                          gamma=5)) # 5 from below parameter tuning
  #####
  # gammas <- 5*(1:20)-4
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=0.5)
  #####
  ### Sample size 2500 (eps=5, p=10) ###
  sample.size <- 2500
  eps <- 1000
  tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                      sample.size=sample.size,
                      gamma=10)) # 10 from below parameter tuning
  #####
  # gammas <- 5*(1:20)-4
  # n.simulations <- 1000
  # evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
  #                                           sample.size, eps,
  #                                           train.proportion=0.5)
}

#####
###### Tuning for eps = Inf
{
### Sample size 200 (eps=Inf, p=10) ###
sample.size <- 200
eps <- Inf
if (neg.weights) {
  gam <- 6
  if (ncol(public.data$X)==4) gam <- 8
}
else gam <- 2
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                            sample.size=sample.size,
                            gamma=gam))
                            # gamma=6)) # 6 for negative weights
                            # gamma=2)) # 2 unnormalized benefit
                            # gamma=0.1)) # 0.1 for trim.B 
                            # gamma=0.1)) # 0.1 normalized benefit

#####
# # gammas <- 1:20/10
# # gammas <- 1:30*5
# # gammas <- 1:20/10 + 1.5
# # gammas <- 1:20/40
# # gammas <- 1:20/100
# gammas <- 1:20
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=2/10)

#####
### Sample size 500 (eps=Inf, p=10) ###
sample.size <- 500
eps <- Inf
if (neg.weights) {
  gam <- 17
  if (ncol(public.data$X)==4) gam <- 15
}
else gam <- 1
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                              sample.size=sample.size,
                              gamma=gam)) 
                              # gamma=17)) # 17 for negative weights
                              # gamma=1)) # 1 unnormalized benefit
                              # gamma=0.2)) # 0.2 for trim.B
                              # gamma=0.2)) # 0.2 normalized benefit

#####
# # gammas <- 1:20/10
# # gammas <- 1:20/40
# gammas <- 1:20
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=0.5)

#####
### Sample size 800 (eps=Inf, p=10) ###
sample.size <- 800
eps <- Inf
if (neg.weights) {
  gam <- 30
  if (ncol(public.data$X)==4) gam <- 55
}
else gam <- 1
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                           sample.size=sample.size,
                           gamma=gam))
                           # gamma=30)) # 30 for negative weights
                           # gamma=1)) # 1 unnormalized benefit
                           # gamma=0.2)) # 0.2 for trim.B
                           # gamma=0.05)) # 0.05 normalized benefit

#####
# # gammas <- 1:30/10
# # gammas <- 1:20/80
# # gammas <- 1:20/40
# # gammas <- 1:30
# gammas <- 2*(5:30)
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=0.5)

#####
### Sample size 1000 (eps=Inf, p=10) ###
sample.size <- 1000
eps <- Inf
if (neg.weights) {
  gam <- 30
  if (ncol(public.data$X)==4) gam <- 50
}
else gam <- 1
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                            sample.size=sample.size,
                            gamma=gam))
                            # gamma=30)) # 30 for negative weights
                            # gamma=1)) # 1 unnormalized benefit
                            # gamma=0.2)) # 0.2 for trim.B
                            # gamma=0.05)) #0.05 normalized benefit
#####
# # gammas <- 1:10/5
# # gammas <- 1:20
# gammas <- 2*(7:40)
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=0.5)
#####
### Sample size 1500 (eps=Inf, p=10) ###
sample.size <- 1500
eps <- Inf
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                              sample.size=sample.size,
                              gamma=2)) # 2 from below parameter tuning
#####
# gammas <- 1:30/10
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=0.5)
#####
### Sample size 2000 (eps=Inf, p=10) ###
sample.size <- 2000
eps <- Inf
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                        sample.size=sample.size,
                        gamma=2)) # 2 from below parameter tuning
#####
# gammas <- 1:30/10
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=0.5)
#####
### Sample size 2500 (eps=Inf, p=10) ###
sample.size <- 2500
eps <- Inf
tuned.reg.constants <- rbind(tuned.reg.constants, data.frame(eps=eps,
                            sample.size=sample.size,
                            gamma=1)) # 1 from below parameter tuning
#####
# gammas <- 1:30/10
# n.simulations <- 1000
# evaluated.gammas <- evaluate.reg.constant(gammas, n.simulations, public.data,
#                                           sample.size, eps,
#                                           train.proportion=0.5)
}

#####
### Best linear model for gamma from eps and sample.size
# lm(log(gamma)~log(eps)+sample.size, data=tuned.reg.constants[1:84,])

#####
### Plots for tuning regularization constant
#####
# ggplot(evaluated.gammas,
#        aes(Regularization.Constant, Mean.Percent.Validation.Accuracy)) +
#   geom_point() +
#   geom_line(linetype="dashed") +
#   geom_errorbar(aes(ymin=Mean.Percent.Validation.Accuracy-qt(0.975, n.simulations-1)*
#                       Std.Dev.Percent.Validation.Accuracy/sqrt(n.simulations),
#                     ymax=Mean.Percent.Validation.Accuracy+qt(0.975, n.simulations-1)*
#                       Std.Dev.Percent.Validation.Accuracy/sqrt(n.simulations)))
# 
# ggplot(evaluated.gammas,
#        aes(Regularization.Constant, Mean.Estimated.Validation.Value)) +
#   geom_point() +
#   geom_line(linetype="dashed") +
#   geom_errorbar(aes(ymin=Mean.Estimated.Validation.Value-qt(0.975, n.simulations-1)*
#                       Std.Dev.Estimated.Validation.Value/sqrt(n.simulations),
#                     ymax=Mean.Estimated.Validation.Value+qt(0.975, n.simulations-1)*
#                       Std.Dev.Estimated.Validation.Value/sqrt(n.simulations)))

###############################################
### Simulation Study for Summary Statistics ###
###############################################
g.factors <- c(.1, .25, .5, 1)
g.factor <- g.factors[4]

simulation.study <- function(n.simulations, N, eps, p, gamma, kernel="linear",
                             D=15, weights.upper.bound=wub, normalize.X=FALSE,
                             n.divisions=1, subset.rate=1){
  # Store list of settings
  settings = list(n.simulations=n.simulations,
                  sample.size=N,
                  epsilon=eps,
                  n.predictors=p,
                  reg.constant=gamma,
                  kernel=kernel,
                  D=D,
                  weights.upper.bound=weights.upper.bound,
                  n.divisions=n.divisions,
                  subset.rate=subset.rate)

  # Run simulations
  prediction.accuracies <- numeric(n.simulations)
  val.func.estimates <- numeric(n.simulations)
  for (i in 1:n.simulations){
    # Generate data
    training.data <- generate.data(N, p, normalize.X=normalize.X)
    validation.data <- generate.data(5000, p, normalize.X=normalize.X)

    # Train model
    wsvm <- train.model(training.data, eps, gamma, kernel, D,
                        weights.upper.bound, n.divisions=n.divisions,
                        subset.rate=subset.rate)

    # Evaluate model
    results <- evaluate.model(wsvm, validation.data)
    prediction.accuracies[i] <- results$percent.accuracy
    val.func.estimates[i] <- results$val.func
  }
  return(list(settings=settings,
              prediction.accuracies=prediction.accuracies,
              val.func.estimates=val.func.estimates))
}

# Simulation settings
n.simulations <- 200
sample.sizes <- c(200, 500, 800, 1000)#, 1500, 2000, 2500)
epsilons <- c(0.1, 0.5, 1, 2, 5, Inf)# 20, 50, 150, 300, 500, 800, 1000, Inf)
n.predictors <- 4 # c(10, 4) # c(10, 20)

# Run simulations
simulation.results <- list()
idx <- 1
for (N in sample.sizes){
  for (eps in epsilons){
    for (p in n.predictors){
      g <- get.reg.constant(tuned.reg.constants, eps,
                            N/n.divisions)/benefit.scale
      if (subset.rate<1){
        g <- get.reg.constant(tuned.reg.constants, eps/(n.divisions*subset.rate),
                              N*subset.rate)/benefit.scale
      }
      # g <- g*(1+g.factor)
      print(paste("Simulating N=", N, ", eps=", eps, ", p=", p,
                  ", gamma=", g, "...", sep=""))
      simulation.results[[idx]] <- simulation.study(n.simulations, N, eps, p,
                                                    gamma=g,
                                                    normalize.X=normal.X,
                                                    n.divisions=n.divisions,
                                                    subset.rate=subset.rate)
      idx <- idx+1
    }
  }
}

#######################
### Analyze Results ###
#######################
# Pull simulation(s) satisfying only desired conditions
get.simulations <- function(simulations, sample.size=NULL, epsilon=NULL,
                            n.predictors=NULL){
  out <- list()
  idx <- 1
  for (i in 1:length(simulations)){
    sim <- simulations[[i]]
    # Check sample size
    if (!is.null(sample.size)){
      if (!(sim$settings$sample.size %in% sample.size)){
        next
      }
    }
    # Check epsilon
    if (!is.null(epsilon)){
      if (!(sim$settings$epsilon %in% epsilon)){
        next
      }
    }
    # Check num predictors
    if (!is.null(n.predictors)){
      if (!(sim$settings$n.predictors %in% n.predictors)){
        next
      }
    }
    # All constraints satisfied
    out[[idx]] <- sim
    idx <- idx+1
  }
  return(out)
}

summarize.simulation <- function(simulation, verbose=FALSE){
  mean.prediction.accuracy <- mean(simulation$prediction.accuracies)
  sd.prediction.accuracy <- sd(simulation$prediction.accuracies)
  mean.val.func.estimate <- mean(simulation$val.func.estimates)
  sd.val.func.estimate <- sd(simulation$val.func.estimates)
  if (verbose){
    print(paste("Setting: N=", simulation$settings$sample.size,
                ", eps=", simulation$settings$epsilon,
                ", p=", simulation$settings$n.predictors,
                ", gamma=", simulation$settings$reg.constant,
                sep=""))
    print("Prediction Accuracy:")
    print(paste("Mean= ", mean.prediction.accuracy,
                ", Sd= ", sd.prediction.accuracy, sep=""))
    print("Estimated Value Function:")
    print(paste("Mean= ", mean.val.func.estimate,
                ", Sd= ", sd.val.func.estimate, sep=""))
  }
  return(list(mean.prediction.accuracy=mean.prediction.accuracy,
              sd.prediction.accuracy=sd.prediction.accuracy,
              mean.val.func.estimate=mean.val.func.estimate,
              sd.val.func.estimate=sd.val.func.estimate))
}

# Use functions to see how results change when fixing sample size and num
#    predictors but varying epsilon
# results.subset <- get.simulations(simulation.results, sample.size=sample.sizes[1],
#                                   n.predictors=n.predictors[1])
# for (i in 1:length(results.subset)){
#   summarize.simulation(results.subset[[i]], TRUE)
#   print("")
# }

# Create data.frame for convenient analysis storage
Ns <- numeric(length(simulation.results))
epss <- numeric(length(simulation.results))
ps <- numeric(length(simulation.results))
gs <- numeric(length(simulation.results))
prediction.accuracy.means <- numeric(length(simulation.results))
prediction.accuracy.sds <- numeric(length(simulation.results))
val.func.estimate.means <- numeric(length(simulation.results))
val.func.estimate.sds <- numeric(length(simulation.results))
for (i in 1:length(simulation.results)){
  sim <- simulation.results[[i]]
  Ns[i] <- sim$settings$sample.size
  epss[i] <- sim$settings$epsilon
  ps[i] <- sim$settings$n.predictors
  gs[i] <- sim$settings$reg.constant
  sim.summary <- summarize.simulation(sim)
  prediction.accuracy.means[i] <- sim.summary$mean.prediction.accuracy
  prediction.accuracy.sds[i] <- sim.summary$sd.prediction.accuracy
  val.func.estimate.means[i] <- sim.summary$mean.val.func.estimate
  val.func.estimate.sds[i] <- sim.summary$sd.val.func.estimate
}
simulation.study.summary.results <- data.frame(
  "Sample Size"=Ns,
  "Epsilon"=epss,
  "Num Predictors"=ps,
  "gamma"=gs,
  "Prediction Accuracy Mean (%)"=prediction.accuracy.means,
  "Prediction Accuracy Std Dev (%)"=prediction.accuracy.sds,
  "Estimated Value Function Mean"=val.func.estimate.means,
  "Estimated Value Function Std Dev"=val.func.estimate.sds,
  check.names = FALSE)

####################
### Plot Results ###
####################
# Convert Inf to number larger than max eps for plotting purposes
simulation.study.summary.results$Epsilon[
  simulation.study.summary.results$Epsilon==Inf] <- 6 # 7000 # 1100 # 1005
simulation.study.summary.results.EPS.INF <- simulation.study.summary.results[
  simulation.study.summary.results$Epsilon==6, # 7000, # 1100, # 1005,
]
simulation.study.summary.results <- simulation.study.summary.results[
  simulation.study.summary.results$Epsilon!=6, # 7000, # 1100, # 1005,
]

# Split on number of predictors p
summary.results.1 <- simulation.study.summary.results[
  simulation.study.summary.results$`Num Predictors`==n.predictors[1],]
summary.results.2 <- simulation.study.summary.results[
  simulation.study.summary.results$`Num Predictors`==n.predictors[2],]
summary.results.1.EPS.INF <- simulation.study.summary.results.EPS.INF[
  simulation.study.summary.results.EPS.INF$`Num Predictors`==n.predictors[1],]
summary.results.2.EPS.INF <- simulation.study.summary.results.EPS.INF[
  simulation.study.summary.results.EPS.INF$`Num Predictors`==n.predictors[2],]

# clipping <- "without clipping"
# if (trim.B) clipping <- "with clipping"
# gam.inc <- as.character(100*g.factor)
# title <- paste("Accuracy (", gam.inc, "% gamma increase ", clipping, ")",
#                sep="")

title <- "Accuracy (without negative weights)"
if (neg.weights) title <- "Accuracy (with negative weights)"

# Plot for first p
ggplot(summary.results.1,
       aes(x=Epsilon, y=`Prediction Accuracy Mean (%)`)) + 
  geom_line(aes(y=`Prediction Accuracy Mean (%)`,
                linetype=factor(`Sample Size`),
                color=factor(`Sample Size`))) +
  geom_point(aes(color=factor(`Sample Size`))) +
  geom_errorbar(aes(ymin=`Prediction Accuracy Mean (%)`-qt(0.975, n.simulations-1)*
                      `Prediction Accuracy Std Dev (%)`/sqrt(n.simulations),
                    ymax=`Prediction Accuracy Mean (%)`+qt(0.975, n.simulations-1)*
                      `Prediction Accuracy Std Dev (%)`/sqrt(n.simulations),
                    color=factor(`Sample Size`)), show.legend=FALSE,
                width=0.15) +
  geom_point(data=summary.results.1.EPS.INF,
             aes(Epsilon, `Prediction Accuracy Mean (%)`,
                 colour=factor(`Sample Size`))) +
  geom_errorbar(data=summary.results.1.EPS.INF,
                aes(ymin=`Prediction Accuracy Mean (%)`-qt(0.975, n.simulations-1)*
                      `Prediction Accuracy Std Dev (%)`/sqrt(n.simulations),
                    ymax=`Prediction Accuracy Mean (%)`+qt(0.975, n.simulations-1)*
                      `Prediction Accuracy Std Dev (%)`/sqrt(n.simulations),
                    color=factor(`Sample Size`)), show.legend=FALSE,
                width=0.15) +
  # scale_x_continuous(trans='log10', breaks=c(10^(-1:3), 7000),
  #                    labels=c(as.character(10^(-1:3)), expression(infinity))) +
  scale_x_continuous(breaks=0:6,
                     labels=c("0", "1", "2", "3", "4", "5", expression(infinity))) +
  guides(color = guide_legend(title = "n"), linetype=guide_legend(title="n")) +
  ylim(40, 100) +
  ylab("Mean Prediction Accuracy (%)") +
  xlab(expression(epsilon)) + 
  theme(text = element_text(size = 16),
        legend.position = c(0.9, 0.2),
        legend.key.width=unit(1,"cm")) +
  ggtitle("Optimal Treatment Prediction Accuracy")
  # scale_x_continuous(breaks=c((5*0:11), c(995, 1000, 1005)),
  #                    labels=c("0", "5", "10", "15", "20", "25", "30", "35",
  #                    "40", "45", "50", "55", "995", "1000", "Infinity")) +
  # scale_x_break(c(55, 995)) +
  # scale_x_continuous(breaks=100*0:11,
  #                    labels=c(as.character(100*0:10), "Infinity")) +
  # ggtitle(paste("Optimal Treatment Prediction Accuracy (p=",
  #               n.predictors[1], ")", sep=""))
  # ggtitle("Accuracy (discretized B)")
  # ggtitle(title)

# # Plot for second p
# ggplot(summary.results.2,
#        aes(Epsilon, `Prediction Accuracy Mean (%)`,
#            colour=factor(`Sample Size`))) + geom_point() +
#   theme(text = element_text(size = 16)) +
#   geom_line(linetype="dashed") +
#   geom_errorbar(aes(ymin=`Prediction Accuracy Mean (%)`-qt(0.975, n.simulations-1)*
#                       `Prediction Accuracy Std Dev (%)`/sqrt(n.simulations),
#                     ymax=`Prediction Accuracy Mean (%)`+qt(0.975, n.simulations-1)*
#                       `Prediction Accuracy Std Dev (%)`/sqrt(n.simulations))) +
#   geom_point(data=summary.results.2.EPS.INF,
#              aes(Epsilon, `Prediction Accuracy Mean (%)`,
#                  colour=factor(`Sample Size`))) +
#   geom_errorbar(data=summary.results.2.EPS.INF,
#                 aes(ymin=`Prediction Accuracy Mean (%)`-qt(0.975, n.simulations-1)*
#                       `Prediction Accuracy Std Dev (%)`/sqrt(n.simulations),
#                     ymax=`Prediction Accuracy Mean (%)`+qt(0.975, n.simulations-1)*
#                       `Prediction Accuracy Std Dev (%)`/sqrt(n.simulations)),
#                 width=0.25) +
#   # scale_x_continuous(trans='log10', breaks=c(10^(-1:3), 7000),
#   #                    labels=c(as.character(10^(-1:3)), "Infinity")) +
#   scale_x_continuous(breaks=0:6,
#                      labels=c("0", "1", "2", "3", "4", "5", "Infinity")) +
#   # scale_x_continuous(breaks=c((5*0:11), c(995, 1000, 1005)),
#   #                    labels=c("0", "5", "10", "15", "20", "25", "30", "35",
#   #                    "40", "45", "50", "55", "995", "1000", "Infinity")) +
#   # scale_x_break(c(55, 995)) +
#   # scale_x_continuous(breaks=100*0:11,
#   #                    labels=c(as.character(100*0:10), "Infinity")) +
#   guides(color = guide_legend(title = "Sample Size")) + ylim(40, 100) +
#   ggtitle(paste("Optimal Treatment Prediction Accuracy (p=",
#                 n.predictors[2], ")", sep=""))
# # ggtitle("Accuracy (discretized B)")
# # ggtitle(title)

title <- "Treatment Value (without negative weights)"
if (neg.weights) title <- "Treatment Value (with negative weights)"

# Plot for first p
ggplot(summary.results.1,
       aes(x=Epsilon, y=`Estimated Value Function Mean`)) +
  geom_line(aes(y=`Estimated Value Function Mean`,
                linetype=factor(`Sample Size`),
                color=factor(`Sample Size`))) +
  geom_point(aes(color=factor(`Sample Size`))) +
  geom_errorbar(aes(ymin=`Estimated Value Function Mean`-qt(0.975, n.simulations-1)*
                      `Estimated Value Function Std Dev`/sqrt(n.simulations),
                    ymax=`Estimated Value Function Mean`+qt(0.975, n.simulations-1)*
                      `Estimated Value Function Std Dev`/sqrt(n.simulations),
                    color=factor(`Sample Size`)), show.legend=FALSE,
                width=0.15) +
  geom_point(data=summary.results.1.EPS.INF,
             aes(Epsilon, `Estimated Value Function Mean`,
                 colour=factor(`Sample Size`))) +
  geom_errorbar(data=summary.results.1.EPS.INF,
                aes(ymin=`Estimated Value Function Mean`-qt(0.975, n.simulations-1)*
                      `Estimated Value Function Std Dev`/sqrt(n.simulations),
                    ymax=`Estimated Value Function Mean`+qt(0.975, n.simulations-1)*
                      `Estimated Value Function Std Dev`/sqrt(n.simulations),
                    color=factor(`Sample Size`)), show.legend=FALSE,
                width=0.15) +
  # scale_x_continuous(trans='log10', breaks=c(10^(-1:3), 7000),
  #                    labels=c(as.character(10^(-1:3)), expression(infinity))) +
  scale_x_continuous(breaks=0:6,
                     labels=c("0", "1", "2", "3", "4", "5", expression(infinity))) +
  guides(color = guide_legend(title = "n"), linetype=guide_legend(title="n")) +
  ylab("Mean Treatment Value") +
  xlab(expression(epsilon)) +
  theme(text = element_text(size = 16),
        legend.position = c(0.9, 0.2),
        legend.key.width=unit(1,"cm")) +
  ggtitle("Empirical Treatment Value")
  # scale_x_continuous(breaks=c((5*0:11), c(995, 1000, 1005)),
  #                    labels=c("0", "5", "10", "15", "20", "25", "30", "35",
  #                    "40", "45", "50", "55", "995", "1000", "Infinity")) +
  # scale_x_break(c(55, 995)) +
  # scale_x_continuous(breaks=100*0:11,
  #                    labels=c(as.character(100*0:10), "Infinity")) +
  # ggtitle("Value (discretized B)")
  # ggtitle(paste("Treatment Value (p=", n.predictors[1], ")", sep=""))
  # ggtitle(title)

# # Plot for second p
# ggplot(summary.results.2,
#        aes(Epsilon, `Estimated Value Function Mean`,
#            colour=factor(`Sample Size`))) + geom_point() +
#   theme(text = element_text(size = 16)) +
#   geom_line(linetype="dashed") +
#   geom_errorbar(aes(ymin=`Estimated Value Function Mean`-qt(0.975, n.simulations-1)*
#                       `Estimated Value Function Std Dev`/sqrt(n.simulations),
#                     ymax=`Estimated Value Function Mean`+qt(0.975, n.simulations-1)*
#                       `Estimated Value Function Std Dev`/sqrt(n.simulations))) +
#   geom_point(data=summary.results.2.EPS.INF,
#              aes(Epsilon, `Estimated Value Function Mean`,
#                  colour=factor(`Sample Size`))) +
#   geom_errorbar(data=summary.results.2.EPS.INF,
#                 aes(ymin=`Estimated Value Function Mean`-qt(0.975, n.simulations-1)*
#                       `Estimated Value Function Std Dev`/sqrt(n.simulations),
#                     ymax=`Estimated Value Function Mean`+qt(0.975, n.simulations-1)*
#                       `Estimated Value Function Std Dev`/sqrt(n.simulations)),
#                 width=0.25) +
#   # scale_x_continuous(trans='log10', breaks=c(10^(-1:3), 7000),
#   #                    labels=c(as.character(10^(-1:3)), "Infinity")) +
#   scale_x_continuous(breaks=0:6,
#                      labels=c("0", "1", "2", "3", "4", "5", "Infinity")) +
#   # scale_x_continuous(breaks=c((5*0:11), c(995, 1000, 1005)),
#   #                    labels=c("0", "5", "10", "15", "20", "25", "30", "35",
#   #                    "40", "45", "50", "55", "995", "1000", "Infinity")) +
#   # scale_x_break(c(55, 995)) +
#   # scale_x_continuous(breaks=100*0:11,
#   #                    labels=c(as.character(100*0:10), "Infinity")) +
#   guides(color = guide_legend(title = "Sample Size")) +
#   # ggtitle("Value (discretized B)")
#   ggtitle(paste("Treatment Value (p=", n.predictors[2], ")", sep=""))


# # Function to get 95% CI for accuracy/value function for give eps and N
# get.CI <- function(eps, N, n.simulations, which.one=c('accuracy', 'value')){
#   # Subset data
#   if (!is.infinite(eps)){
#     results <- simulation.study.summary.results[
#       simulation.study.summary.results$`Sample Size`==N,]
#     results <- results[results$Epsilon==eps,]
#   } else{
#     results <- simulation.study.summary.results.EPS.INF[
#       simulation.study.summary.results.EPS.INF$`Sample Size`==N,]
#   }
#   
#   if (which.one=='accuracy'){
#     mean.info <- results$`Prediction Accuracy Mean (%)`
#     std.info <- results$`Prediction Accuracy Std Dev (%)`
#   } else if (which.one=='value'){
#     mean.info <- results$`Estimated Value Function Mean`
#     std.info <- results$`Estimated Value Function Std Dev`
#   }
#   
#   out <- paste("(", round(mean.info - qt(0.975, n.simulations-1)*
#                  std.info/sqrt(n.simulations), 2), ", ",
#                round(mean.info + qt(0.975, n.simulations-1)*
#                  std.info/sqrt(n.simulations), 2), ")", sep="")
#   out
# }
# 
# acc.CIs <- c()
# for (eps in c(0.1, 0.5, 1, 2, 5, Inf)){
#   for (N in c(200, 500, 800, 1000)){
#     acc.CIs <- c(acc.CIs, get.CI(eps, N, n.simulations, 'accuracy'))
#   }
# }
# 
# val.CIs <- c()
# for (eps in c(0.1, 0.5, 1, 2, 5, Inf)){
#   for (N in c(200, 500, 800, 1000)){
#     val.CIs <- c(val.CIs, get.CI(eps, N, n.simulations, 'value'))
#   }
# }
