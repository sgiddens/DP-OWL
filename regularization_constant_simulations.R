library(DPpack)
library(ggplot2)
library(truncnorm)
library(dplyr)
library(latex2exp)

#####################
### Generate data ###
#####################
generate.data <- function(N, p, theta, propensities=c(0.5, 0.5),
                          normalize.X=FALSE, X.dist='uniform',
                          norm.mus=0, norm.sds=1){
  # Generate X
  X <- matrix(NaN, nrow=N, ncol=p)
  for (i in 1:N){
    if (X.dist=='uniform'){
      tmp <- runif(p)
    } else if (X.dist=='truncnorm'){
      tmp <- rtruncnorm(p, a=0, b=1, mean=norm.mus, sd=norm.sds)
    }
    if (normalize.X) tmp <- tmp/sqrt(p+1) # Normalize X
    X[i, ] <- tmp
  }
  X <- data.frame(X)

  # Generate A
  A <- sample(c(-1, 1), size=N, replace=TRUE, prob=propensities)

  # Ground truth function
  f <- function(x) cbind(1, as.matrix(x))%*%as.matrix(theta[1:(ncol(x)+1)])

  # Generate B
  mus <- 0.01 + 0.02*X[,4] + 3*A*f(X)
  std.dev <- 0.5
  B <- numeric(N)
  for (i in 1:N){
    B[i] <- rnorm(1, mean=mus[i], sd=std.dev)
  }
  if (min(B) < 0) B <- B + abs(min(B)) + 0.001

  return(list(X=X, A=A, B=B, f=f))
}

# Private data properties
p <- 10
theta <- c(1, 1, -1.8, -2.2)
theta <- c(theta, rep(0, p-length(theta)))

################################
### Train wSVM model with DP ###
################################
train.model <- function(training.data, eps, gamma, kernel, D,
                        weights.upper.bound, propensity=0.5, add.bias=TRUE,
                        normalize.X=FALSE){
  # Redefine target variable to align with svmDP requirements
  y <- training.data$A
  y[y<0] <- 0
  # Bounds on training data (0 and 1 since drawn from U(0, 1))
  p <- ncol(training.data$X)
  upper.bounds <- rep(1, p)
  if (normalize.X) upper.bounds <- upper.bounds/sqrt(p+1)
  lower.bounds <- rep(0, p)
  # Weights from B and propensity score
  weights <- training.data$B/propensity

  # set.seed(NULL)
  # Construct object
  wsvm <- svmDP$new("l2", eps, gamma, perturbation.method="output", kernel,
                    D=D)

  # Train model
  wsvm$fit(training.data$X, y, upper.bounds, lower.bounds,
           add.bias, weights, weights.upper.bound)
  return(wsvm)
}

###########################
### Evaluate wSVM model ###
###########################
# Estimate value function
estimate.value.function <- function(trained.model, val.X, val.A, val.B,
                                    propensity, add.bias){
  val.A.hat <- trained.model$predict(val.X, add.bias)
  val.A.hat[val.A.hat==0] <- -1
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

################################################
### Plots for tuning regularization constant ###
################################################
plot.gamma.acc <- function(evaluated.gammas){
  ggplot(evaluated.gammas,
         aes(Regularization.Constant, Mean.Percent.Validation.Accuracy)) +
    geom_point() +
    geom_line(linetype="dashed") +
    geom_errorbar(aes(ymin=Mean.Percent.Validation.Accuracy-qt(0.975, n.simulations-1)*
                        Std.Dev.Percent.Validation.Accuracy/sqrt(n.simulations),
                      ymax=Mean.Percent.Validation.Accuracy+qt(0.975, n.simulations-1)*
                        Std.Dev.Percent.Validation.Accuracy/sqrt(n.simulations)))
}

plot.gamma.val.func <- function(evaluated.gammas){
  ggplot(evaluated.gammas,
         aes(Regularization.Constant, Mean.Estimated.Validation.Value)) +
    geom_point() +
    geom_line(linetype="dashed") +
    geom_errorbar(aes(ymin=Mean.Estimated.Validation.Value-qt(0.975, n.simulations-1)*
                        Std.Dev.Estimated.Validation.Value/sqrt(n.simulations),
                      ymax=Mean.Estimated.Validation.Value+qt(0.975, n.simulations-1)*
                        Std.Dev.Estimated.Validation.Value/sqrt(n.simulations)))
}

############################################
##### Simulations for comparison plots #####
############################################
get.comparison.data <- function(epsilons, sample.sizes, public.data, wub,
                                gamma.grid=20*(1:21)-19, n.simulations=1000,
                                seed=42){
  set.seed(seed)
  results <- data.frame()
  for (eps in epsilons){
    for (sample.size in sample.sizes){
      print("----------------------------")
      print(paste0("epsilon=", eps, "; sample size=", sample.size))
      tmp <- evaluate.reg.constant(gamma.grid, n.simulations,
                                   public.data, sample.size, eps,
                                   weights.upper.bound=wub)
      tmp$epsilon <- eps
      tmp$sample.size <- sample.size
      results <- rbind(results, tmp)
    }
  }
  return(results)
}

############################
##### Comparison plots #####
############################
plot.comparison.data <- function(epsilons, sample.sizes, comparison.data,
                                 n.simulations=1000, subfolder=""){
  title.size <- 50
  axis.title.size <- 45
  axis.text.size <- 45
  axis.tick.size <- 5
  for (i in 1:length(epsilons)){
    for (j in 1:length(sample.sizes)){
      eps <- epsilons[i]
      ss <- sample.sizes[j]
      evaluated.gammas <- comparison.data %>% filter(epsilon==eps,
                                                     sample.size==ss)

      ggp <- ggplot(evaluated.gammas, aes(Regularization.Constant,
                                          Mean.Percent.Validation.Accuracy)) +
        geom_point() +
        geom_line(linetype='dashed') +
        geom_errorbar(
          aes(ymin=Mean.Percent.Validation.Accuracy-qnorm(0.975)*
                Std.Dev.Percent.Validation.Accuracy/sqrt(n.simulations),
              ymax=Mean.Percent.Validation.Accuracy+qnorm(0.975)*
                Std.Dev.Percent.Validation.Accuracy/sqrt(n.simulations)
          )
        )
      ggp <- ggp +
        theme(axis.ticks=element_line(linewidth=axis.tick.size)) +
        scale_y_continuous(limits=c(48.5,101.5), expand=c(0,0))


      height <- 7
      width <- 7
      if (i==1 & j==1){
        ggp <- ggp +
          theme(text = element_text(size = axis.text.size),
                legend.position = c(0.5, 0.7),
                legend.key.width=unit(1,"cm"))
      } else{
        ggp <- ggp + theme(legend.position="none")
      }
      if (i==1){ # Add title to top row
        ggp <- ggp + ggtitle(TeX(paste0("$n=", ss, "$"))) +
          theme(plot.title=element_text(size=title.size, hjust=0.5))
      }
      if (j==1){ # Add y-axis label to left column
        ggp <- ggp + ylab(TeX(paste0("$\\epsilon = ", eps, "$"))) +
          theme(axis.title.y=element_text(size=axis.title.size),
                axis.text.y=element_text(size=axis.text.size))
        width <- 8.6 # Accounts for axis labels
      } else{
        ggp <- ggp + theme(axis.title.y=element_blank(),
                           axis.text.y=element_blank())
      }
      if (i==length(epsilons)){ # Add x-axis label to bottom row
        ggp <- ggp + xlab(TeX("$\\gamma$")) +
          theme(axis.title.x = element_text(size=axis.title.size),
                axis.text.x = element_text(size=axis.text.size))
        height <- 8.1 # Accounts for axis labels
      } else{
        ggp <- ggp + theme(axis.title.x = element_blank(),
                           axis.text.x = element_blank())
      }

      str_eps <- gsub("\\.", "-", as.character(eps))
      ggsave(ggp, file=paste0("RegConstFigures/", subfolder,
                              "eps_", str_eps, "_N_", ss, ".pdf"),
             height=height, width=width)
    }
  }
}

plot.comparison.data.N <- function(epsilons, sample.sizes, comparison.data,
                                 n.simulations=1000, plot.vert.lines=FALSE,
                                 subfolder=""){
  title.size <- 50
  axis.title.size <- 45
  axis.text.size <- 45
  axis.tick.size <- 5
  for (i in 1:length(epsilons)){
    for (j in 1:length(sample.sizes)){
      eps <- epsilons[i]
      ss <- sample.sizes[j]
      evaluated.gammas <- comparison.data %>% filter(epsilon==eps,
                                                     sample.size==ss)

      ggp <- ggplot(evaluated.gammas, aes(Regularization.Constant,
                                          Mean.Percent.Validation.Accuracy),
                    color=factor(N)
                    ) +
        geom_point(
          aes(color=factor(N))
          ) +
        geom_line(
          aes(color=factor(N)),
          linetype='dashed') +
        geom_errorbar(
          aes(ymin=Mean.Percent.Validation.Accuracy-qnorm(0.975)*
                Std.Dev.Percent.Validation.Accuracy/sqrt(n.simulations),
              ymax=Mean.Percent.Validation.Accuracy+qnorm(0.975)*
                Std.Dev.Percent.Validation.Accuracy/sqrt(n.simulations),
              color=factor(N)
              )
        ) +
        scale_color_manual(values=c("#7CAE00", "#00BFC4", "black", "#C77CFF"),
                           labels=c("100", "500", "1000 (baseline)", "2500"))

      ggp <- ggp +
        theme(axis.ticks=element_line(linewidth=axis.tick.size)) +
        scale_y_continuous(limits=c(48.5,101.5), expand=c(0,0))


      height <- 7
      width <- 7
      if (i==1 & j==1){
        ggp <- ggp + guides(color = guide_legend(title = TeX("$n_0$"))) +
          theme(text = element_text(size = axis.text.size),
                legend.position = c(0.5, 0.7),
                legend.key.width=unit(1,"cm"))
      } else{
        ggp <- ggp + theme(legend.position="none")
      }
      if (i==1){ # Add title to top row
        ggp <- ggp + ggtitle(TeX(paste0("$n=", ss, "$"))) +
          theme(plot.title=element_text(size=title.size, hjust=0.5))
      }
      if (j==1){ # Add y-axis label to left column
        ggp <- ggp + ylab(TeX(paste0("$\\epsilon = ", eps, "$"))) +
          theme(axis.title.y=element_text(size=axis.title.size),
                axis.text.y=element_text(size=axis.text.size))
        width <- 8.6 # Accounts for axis labels
      } else{
        ggp <- ggp + theme(axis.title.y=element_blank(),
                           axis.text.y=element_blank())
      }
      if (i==length(epsilons)){ # Add x-axis label to bottom row
        ggp <- ggp + xlab(TeX("$\\gamma$")) +
          theme(axis.title.x = element_text(size=axis.title.size),
                axis.text.x = element_text(size=axis.text.size))
        height <- 8.1 # Accounts for axis labels
      } else{
        ggp <- ggp + theme(axis.title.x = element_blank(),
                           axis.text.x = element_blank())
      }

      left.x <- 1
      right.x <- 601
      if (plot.vert.lines){
        evaluated.gammas <- evaluated.gammas %>% filter(N==1000)
        fun <- approxfun(evaluated.gammas$Regularization.Constant,
                         evaluated.gammas$Mean.Percent.Validation.Accuracy)
        optimum <- as.vector(evaluated.gammas[
          which.max(evaluated.gammas$Mean.Percent.Validation.Accuracy), c(1,2)])
        max.minus.5 <- optimum$Mean.Percent.Validation.Accuracy - 5
        acc.1 <- evaluated.gammas[
          evaluated.gammas$Regularization.Constant==1,]$Mean.Percent.Validation.Accuracy
        acc.601 <- evaluated.gammas[
          evaluated.gammas$Regularization.Constant==601,]$Mean.Percent.Validation.Accuracy
        if (acc.1 < max.minus.5){
          search.range <- 1:optimum$Regularization.Constant
          obj.fun <- function(gamma) abs(fun(gamma) - max.minus.5)
          opt.idx <- which.min(obj.fun(search.range))
          left.x <- search.range[opt.idx]
        }
        if (acc.601 < max.minus.5){
          search.range <- optimum$Regularization.Constant:601
          obj.fun <- function(gamma) abs(fun(gamma) - max.minus.5)
          opt.idx <- which.min(obj.fun(search.range))
          right.x <- search.range[opt.idx]
        }
        left.y <- fun(left.x)
        right.y <- fun(right.x)

        ggp <- ggp +
          geom_segment(aes(x = left.x, y = 48.5, xend = left.x, yend = left.y),
                       color="red") +
          geom_segment(aes(x = right.x, y = 48.5, xend = right.x, yend = right.y),
                       color="red")
      }

      str_eps <- gsub("\\.", "-", as.character(eps))
      ggsave(ggp, file=paste0("RegConstFigures/", subfolder,
                              "eps_", str_eps, "_N_", ss, ".pdf"),
             height=height, width=width)
    }
  }
}

plot.comparison.data.diff.theta <- function(epsilons, sample.sizes,
                                           comparison.data, n.simulations=1000,
                                           plot.vert.lines=FALSE, subfolder=""){
  title.size <- 50
  axis.title.size <- 45
  axis.text.size <- 45
  axis.tick.size <- 5
  for (i in 1:length(epsilons)){
    for (j in 1:length(sample.sizes)){
      eps <- epsilons[i]
      ss <- sample.sizes[j]
      evaluated.gammas <- comparison.data %>% filter(epsilon==eps,
                                                     sample.size==ss)

      ggp <- ggplot(evaluated.gammas, aes(Regularization.Constant,
                                          Mean.Percent.Validation.Accuracy),
                    color=factor(theta)
      ) +
        geom_point(
          aes(color=factor(theta))
        ) +
        geom_line(
          aes(color=factor(theta)),
          linetype='dashed') +
        geom_errorbar(
          aes(ymin=Mean.Percent.Validation.Accuracy-qnorm(0.975)*
                Std.Dev.Percent.Validation.Accuracy/sqrt(n.simulations),
              ymax=Mean.Percent.Validation.Accuracy+qnorm(0.975)*
                Std.Dev.Percent.Validation.Accuracy/sqrt(n.simulations),
              color=factor(theta)
          )
        ) +
        scale_color_manual(values=c("black", "#7CAE00", "#C77CFF"),
                           labels=c(TeX("$\\theta_0$ (baseline)"),
                                    TeX("$\\theta_1$"), TeX("$\\theta_2$")))
      ggp <- ggp +
        theme(axis.ticks=element_line(linewidth=axis.tick.size)) +
        scale_y_continuous(limits=c(48.5, 101.5), expand=c(0, 0))


      height <- 7
      width <- 7
      if (i==1 & j==1){
        ggp <- ggp + guides(color = guide_legend(title = TeX("$\\theta$"))) +
          theme(text = element_text(size = axis.text.size),
                legend.position = c(0.5, 0.7),
                legend.key.width=unit(1,"cm"))
      } else{
        ggp <- ggp + theme(legend.position="none")
      }
      if (i==1){ # Add title to top row
        ggp <- ggp + ggtitle(TeX(paste0("$n=", ss, "$"))) +
          theme(plot.title=element_text(size=title.size, hjust=0.5))
      }
      if (j==1){ # Add y-axis label to left column
        ggp <- ggp + ylab(TeX(paste0("$\\epsilon = ", eps, "$"))) +
          theme(axis.title.y=element_text(size=axis.title.size),
                axis.text.y=element_text(size=axis.text.size))
        width <- 8.6 # Accounts for axis labels
      } else{
        ggp <- ggp + theme(axis.title.y=element_blank(),
                           axis.text.y=element_blank())
      }
      if (i==length(epsilons)){ # Add x-axis label to bottom row
        ggp <- ggp + xlab(TeX("$\\gamma$")) +
          theme(axis.title.x = element_text(size=axis.title.size),
                axis.text.x = element_text(size=axis.text.size))
        height <- 8.1 # Accounts for axis labels
      } else{
        ggp <- ggp + theme(axis.title.x = element_blank(),
                           axis.text.x = element_blank())
      }

      left.x <- 1
      right.x <- 601
      if (plot.vert.lines){
        evaluated.gammas <- evaluated.gammas %>% filter(theta=="theta0")
        fun <- approxfun(evaluated.gammas$Regularization.Constant,
                         evaluated.gammas$Mean.Percent.Validation.Accuracy)
        optimum <- as.vector(evaluated.gammas[
          which.max(evaluated.gammas$Mean.Percent.Validation.Accuracy), c(1,2)])
        max.minus.5 <- optimum$Mean.Percent.Validation.Accuracy - 5
        acc.1 <- evaluated.gammas[
          evaluated.gammas$Regularization.Constant==1,]$Mean.Percent.Validation.Accuracy
        acc.601 <- evaluated.gammas[
          evaluated.gammas$Regularization.Constant==601,]$Mean.Percent.Validation.Accuracy
        if (acc.1 < max.minus.5){
          search.range <- 1:optimum$Regularization.Constant
          obj.fun <- function(gamma) abs(fun(gamma) - max.minus.5)
          opt.idx <- which.min(obj.fun(search.range))
          left.x <- search.range[opt.idx]
        }
        if (acc.601 < max.minus.5){
          search.range <- optimum$Regularization.Constant:601
          obj.fun <- function(gamma) abs(fun(gamma) - max.minus.5)
          opt.idx <- which.min(obj.fun(search.range))
          right.x <- search.range[opt.idx]
        }
        left.y <- fun(left.x)
        right.y <- fun(right.x)

        ggp <- ggp +
          geom_segment(aes(x = left.x, y = 48.5, xend = left.x, yend = left.y),
                       color="red") +
          geom_segment(aes(x = right.x, y = 48.5, xend = right.x, yend = right.y),
                       color="red")
      }

      str_eps <- gsub("\\.", "-", as.character(eps))
      ggsave(ggp, file=paste0("RegConstFigures/", subfolder,
                              "eps_", str_eps, "_N_", ss, ".pdf"),
             height=height, width=width)
    }
  }
}

plot.comparison.data.diff.dist <- function(epsilons, sample.sizes,
                                           comparison.data, n.simulations=1000,
                                           plot.vert.lines=FALSE, subfolder=""){
  title.size <- 50
  axis.title.size <- 45
  axis.text.size <- 45
  axis.tick.size <- 5
  for (i in 1:length(epsilons)){
    for (j in 1:length(sample.sizes)){
      eps <- epsilons[i]
      ss <- sample.sizes[j]
      evaluated.gammas <- comparison.data %>% filter(epsilon==eps,
                                                     sample.size==ss)

      ggp <- ggplot(evaluated.gammas, aes(Regularization.Constant,
                                          Mean.Percent.Validation.Accuracy),
                    color=factor(Distribution)
      ) +
        geom_point(
          aes(color=factor(Distribution))
        ) +
        geom_line(
          aes(color=factor(Distribution)),
          linetype='dashed') +
        geom_errorbar(
          aes(ymin=Mean.Percent.Validation.Accuracy-qnorm(0.975)*
                Std.Dev.Percent.Validation.Accuracy/sqrt(n.simulations),
              ymax=Mean.Percent.Validation.Accuracy+qnorm(0.975)*
                Std.Dev.Percent.Validation.Accuracy/sqrt(n.simulations),
              color=factor(Distribution)
          )
        ) +
        scale_color_manual(values=c("#7CAE00", "black"),
                           labels=c(TeX("$X\\sim TN(\\mu, \\sigma^2, 0, 1)$"),
                                    TeX("$X\\sim$Uniform$(0,1)$")))
      ggp <- ggp +
        theme(axis.ticks=element_line(linewidth=axis.tick.size)) +
        scale_y_continuous(limits=c(48.5, 101.5), expand=c(0, 0))


      height <- 7
      width <- 7
      if (i==1 & j==1){
        ggp <- ggp + guides(color = guide_legend(title = "Distribution")) +
          theme(text = element_text(size = axis.text.size),
                legend.position = c(0.5, 0.7),
                legend.key.width=unit(1,"cm"))
      } else{
        ggp <- ggp + theme(legend.position="none")
      }
      if (i==1){ # Add title to top row
        ggp <- ggp + ggtitle(TeX(paste0("$n=", ss, "$"))) +
          theme(plot.title=element_text(size=title.size, hjust=0.5))
      }
      if (j==1){ # Add y-axis label to left column
        ggp <- ggp + ylab(TeX(paste0("$\\epsilon = ", eps, "$"))) +
          theme(axis.title.y=element_text(size=axis.title.size),
                axis.text.y=element_text(size=axis.text.size))
        width <- 8.6 # Accounts for axis labels
      } else{
        ggp <- ggp + theme(axis.title.y=element_blank(),
                           axis.text.y=element_blank())
      }
      if (i==length(epsilons)){ # Add x-axis label to bottom row
        ggp <- ggp + xlab(TeX("$\\gamma$")) +
          theme(axis.title.x = element_text(size=axis.title.size),
                axis.text.x = element_text(size=axis.text.size))
        height <- 8.1 # Accounts for axis labels
      } else{
        ggp <- ggp + theme(axis.title.x = element_blank(),
                           axis.text.x = element_blank())
      }

      left.x <- 1
      right.x <- 601
      if (plot.vert.lines){
        evaluated.gammas <- evaluated.gammas %>% filter(Distribution=="Uniform")
        fun <- approxfun(evaluated.gammas$Regularization.Constant,
                         evaluated.gammas$Mean.Percent.Validation.Accuracy)
        optimum <- as.vector(evaluated.gammas[
          which.max(evaluated.gammas$Mean.Percent.Validation.Accuracy), c(1,2)])
        max.minus.5 <- optimum$Mean.Percent.Validation.Accuracy - 5
        acc.1 <- evaluated.gammas[
          evaluated.gammas$Regularization.Constant==1,]$Mean.Percent.Validation.Accuracy
        acc.601 <- evaluated.gammas[
          evaluated.gammas$Regularization.Constant==601,]$Mean.Percent.Validation.Accuracy
        if (acc.1 < max.minus.5){
          search.range <- 1:optimum$Regularization.Constant
          obj.fun <- function(gamma) abs(fun(gamma) - max.minus.5)
          opt.idx <- which.min(obj.fun(search.range))
          left.x <- search.range[opt.idx]
        }
        if (acc.601 < max.minus.5){
          search.range <- optimum$Regularization.Constant:601
          obj.fun <- function(gamma) abs(fun(gamma) - max.minus.5)
          opt.idx <- which.min(obj.fun(search.range))
          right.x <- search.range[opt.idx]
        }
        left.y <- fun(left.x)
        right.y <- fun(right.x)

        ggp <- ggp +
          geom_segment(aes(x = left.x, y = 48.5, xend = left.x, yend = left.y),
                       color="red") +
          geom_segment(aes(x = right.x, y = 48.5, xend = right.x, yend = right.y),
                       color="red")
      }

      str_eps <- gsub("\\.", "-", as.character(eps))
      ggsave(ggp, file=paste0("RegConstFigures/", subfolder,
                              "eps_", str_eps, "_N_", ss, ".pdf"),
             height=height, width=width)
    }
  }
}

plot.comparison.data.P <- function(epsilons, sample.sizes, comparison.data,
                                   n.simulations=1000, plot.vert.lines=FALSE,
                                   subfolder=""){
  title.size <- 50
  axis.title.size <- 45
  axis.text.size <- 45
  axis.tick.size <- 5
  for (i in 1:length(epsilons)){
    for (j in 1:length(sample.sizes)){
      eps <- epsilons[i]
      ss <- sample.sizes[j]
      evaluated.gammas <- comparison.data %>% filter(epsilon==eps,
                                                     sample.size==ss)

      ggp <- ggplot(evaluated.gammas, aes(Regularization.Constant,
                                          Mean.Percent.Validation.Accuracy),
                    color=factor(p)
      ) +
        geom_point(
          aes(color=factor(p))
        ) +
        geom_line(
          aes(color=factor(p)),
          linetype='dashed') +
        geom_errorbar(
          aes(ymin=Mean.Percent.Validation.Accuracy-qnorm(0.975)*
                Std.Dev.Percent.Validation.Accuracy/sqrt(n.simulations),
              ymax=Mean.Percent.Validation.Accuracy+qnorm(0.975)*
                Std.Dev.Percent.Validation.Accuracy/sqrt(n.simulations),
              color=factor(p)
          )
        ) +
        scale_color_manual(values=c("#7CAE00", "black", "#00BFC4", "#C77CFF"),
                           labels=c("2", "4 (baseline)", "6", "8"))

      ggp <- ggp +
        theme(axis.ticks=element_line(linewidth=axis.tick.size)) +
        scale_y_continuous(limits=c(48.5,101.5), expand=c(0,0))


      height <- 7
      width <- 7
      if (i==1 & j==1){
        ggp <- ggp + guides(color = guide_legend(title = "p")) +
          theme(text = element_text(size = axis.text.size),
                legend.position = c(0.5, 0.7),
                legend.key.width=unit(1,"cm"))
      } else{
        ggp <- ggp + theme(legend.position="none")
      }
      if (i==1){ # Add title to top row
        ggp <- ggp + ggtitle(TeX(paste0("$n=", ss, "$"))) +
          theme(plot.title=element_text(size=title.size, hjust=0.5))
      }
      if (j==1){ # Add y-axis label to left column
        ggp <- ggp + ylab(TeX(paste0("$\\epsilon = ", eps, "$"))) +
          theme(axis.title.y=element_text(size=axis.title.size),
                axis.text.y=element_text(size=axis.text.size))
        width <- 8.6 # Accounts for axis labels
      } else{
        ggp <- ggp + theme(axis.title.y=element_blank(),
                           axis.text.y=element_blank())
      }
      if (i==length(epsilons)){ # Add x-axis label to bottom row
        ggp <- ggp + xlab(TeX("$\\gamma$")) +
          theme(axis.title.x = element_text(size=axis.title.size),
                axis.text.x = element_text(size=axis.text.size))
        height <- 8.1 # Accounts for axis labels
      } else{
        ggp <- ggp + theme(axis.title.x = element_blank(),
                           axis.text.x = element_blank())
      }

      left.x <- 1
      right.x <- 601
      if (plot.vert.lines){
        evaluated.gammas <- evaluated.gammas %>% filter(p==4)
        fun <- approxfun(evaluated.gammas$Regularization.Constant,
                         evaluated.gammas$Mean.Percent.Validation.Accuracy)
        optimum <- as.vector(evaluated.gammas[
          which.max(evaluated.gammas$Mean.Percent.Validation.Accuracy), c(1,2)])
        max.minus.5 <- optimum$Mean.Percent.Validation.Accuracy - 5
        acc.1 <- evaluated.gammas[
          evaluated.gammas$Regularization.Constant==1,]$Mean.Percent.Validation.Accuracy
        acc.601 <- evaluated.gammas[
          evaluated.gammas$Regularization.Constant==601,]$Mean.Percent.Validation.Accuracy
        if (acc.1 < max.minus.5){
          search.range <- 1:optimum$Regularization.Constant
          obj.fun <- function(gamma) abs(fun(gamma) - max.minus.5)
          opt.idx <- which.min(obj.fun(search.range))
          left.x <- search.range[opt.idx]
        }
        if (acc.601 < max.minus.5){
          search.range <- optimum$Regularization.Constant:601
          obj.fun <- function(gamma) abs(fun(gamma) - max.minus.5)
          opt.idx <- which.min(obj.fun(search.range))
          right.x <- search.range[opt.idx]
        }
        left.y <- fun(left.x)
        right.y <- fun(right.x)

        ggp <- ggp +
          geom_segment(aes(x = left.x, y = 48.5, xend = left.x, yend = left.y),
                       color="red") +
          geom_segment(aes(x = right.x, y = 48.5, xend = right.x, yend = right.y),
                       color="red")
      }

      str_eps <- gsub("\\.", "-", as.character(eps))
      ggsave(ggp, file=paste0("RegConstFigures/", subfolder,
                              "eps_", str_eps, "_N_", ss, ".pdf"),
             height=height, width=width)
    }
  }
}

############################
### CAMERA-READY FIGURES ###
############################
### NOTE: Before running the below simulations, you will need to create folders
### "RegConstData/" and "RegConstFigures/" to hold the results. Appropriate 
### subfolders within RegConstFigures/ should also be added to hold the 
### simulation figure results. Subfolder names should match "subfolder"argument 
### in each plotting function.
epsilons <- c(0.5, 1, 5, 20, 150)
sample.sizes <- c(200, 500, 1000, 2000)

# Different n_0
comparison.data1 <- read.csv("RegConstData/match_private_N100.csv")
comparison.data2 <- read.csv("RegConstData/match_private_N500.csv")
comparison.data3 <- read.csv("RegConstData/match_private.csv")
comparison.data4 <- read.csv("RegConstData/match_private_N2500.csv")
comparison.data <- rbind(comparison.data1 %>% mutate(N=100),
                         comparison.data2 %>% mutate(N=500),
                         comparison.data3 %>% mutate(N=1000),
                         comparison.data4 %>% mutate(N=2500))
plot.comparison.data.N(epsilons, sample.sizes, comparison.data,
                       plot.vert.lines=TRUE,
                       subfolder="SupplMaterialsN/")

# Different theta
comparison.data1 <- read.csv("RegConstData/match_private.csv")
comparison.data2 <- read.csv("RegConstData/different_theta.csv")
comparison.data3 <- read.csv("RegConstData/larger_theta.csv")
comparison.data <- rbind(comparison.data1 %>% mutate(theta="theta0"),
                         comparison.data2 %>% mutate(theta="theta1"),
                         comparison.data3 %>% mutate(theta="theta2"))
plot.comparison.data.diff.theta(epsilons, sample.sizes, comparison.data,
                               plot.vert.lines=TRUE,
                               subfolder="SupplMaterialsTheta/")

# Different distribution
comparison.data1 <- read.csv("RegConstData/match_private.csv")
comparison.data2 <- read.csv("RegConstData/truncnorm.csv")
comparison.data <- rbind(comparison.data1 %>% mutate(Distribution="Uniform"),
                         comparison.data2 %>% mutate(Distribution="Normal"))
plot.comparison.data.diff.dist(epsilons, sample.sizes, comparison.data,
                               plot.vert.lines=TRUE,
                               subfolder="SupplMaterialsDist/")

# Different non-zero coefficients (p)
comparison.data1 <- read.csv("RegConstData/different_p2.csv")
comparison.data2 <- read.csv("RegConstData/match_private.csv")
comparison.data3 <- read.csv("RegConstData/different_p6.csv")
comparison.data4 <- read.csv("RegConstData/different_p8.csv")
comparison.data <- rbind(comparison.data1 %>% mutate(p=2),
                         comparison.data2 %>% mutate(p=4),
                         comparison.data3 %>% mutate(p=6),
                         comparison.data4 %>% mutate(p=8))
plot.comparison.data.P(epsilons, sample.sizes, comparison.data,
                       plot.vert.lines=TRUE,
                       subfolder="SupplMaterialsP/")

###################
### SIMULATIONS ###
###################

######## PUBLIC DATA THAT MATCHES PRIVATE DATA #######
set.seed(42)
p <- 10
theta <- c(1, 1, 1, -1.8, -2.2)
theta <- c(theta, rep(0, p+1-length(theta)))
public.data <- generate.data(1000, p, theta)
public.data$X <- public.data$X[,1:4]

gamma.grid <- 20*(1:31)-19
wub <- 30

# comparison.data <- get.comparison.data(epsilons, sample.sizes, public.data, wub,
#                                        gamma.grid)
# write.csv(comparison.data, "RegConstData/match_private.csv", row.names=FALSE)
comparison.data <- read.csv("RegConstData/match_private.csv")
plot.comparison.data(epsilons, sample.sizes, comparison.data,
                     subfolder="MatchPrivate/")

######## PUBLIC DATA OF SMALLER SIZE (N=500) #######
set.seed(42)
p <- 10
theta <- c(1, 1, 1, -1.8, -2.2)
theta <- c(theta, rep(0, p+1-length(theta)))
public.data <- generate.data(500, p, theta)
public.data$X <- public.data$X[,1:4]

gamma.grid <- 20*(1:31)-19
wub <- 30

# comparison.data <- get.comparison.data(epsilons, sample.sizes, public.data, wub,
#                                        gamma.grid)
# write.csv(comparison.data, "RegConstData/match_private_N500.csv", row.names=FALSE)
comparison.data <- read.csv("RegConstData/match_private_N500.csv")
plot.comparison.data(epsilons, sample.sizes, comparison.data,
                     subfolder="MatchPrivateN500/")

######## PUBLIC DATA OF LARGER SIZE (N=2500) #######
set.seed(42)
p <- 10
theta <- c(1, 1, 1, -1.8, -2.2)
theta <- c(theta, rep(0, p+1-length(theta)))
public.data <- generate.data(2500, p, theta)
public.data$X <- public.data$X[,1:4]

gamma.grid <- 20*(1:31)-19
wub <- 30

# comparison.data <- get.comparison.data(epsilons, sample.sizes, public.data, wub,
#                                        gamma.grid)
# write.csv(comparison.data, "RegConstData/match_private_N2500.csv", row.names=FALSE)
comparison.data <- read.csv("RegConstData/match_private_N2500.csv")
plot.comparison.data(epsilons, sample.sizes, comparison.data,
                     subfolder="MatchPrivateN2500/")

######## PUBLIC DATA VERY SMALL SIZE (N=100) #######
set.seed(42)
p <- 10
theta <- c(1, 1, 1, -1.8, -2.2)
theta <- c(theta, rep(0, p+1-length(theta)))
public.data <- generate.data(100, p, theta)
public.data$X <- public.data$X[,1:4]

gamma.grid <- 20*(1:31)-19
wub <- 30

# comparison.data <- get.comparison.data(epsilons, sample.sizes, public.data, wub,
#                                        gamma.grid)
# write.csv(comparison.data, "RegConstData/match_private_N100.csv", row.names=FALSE)
comparison.data <- read.csv("RegConstData/match_private_N100.csv")
plot.comparison.data(epsilons, sample.sizes, comparison.data,
                     subfolder="MatchPrivateN100/")

######## PUBLIC DATA WITH DIFFERENT THETA #######
set.seed(42)
p <- 10
theta <- c(1, -1, -1, 1.5, -1.5)
theta <- c(theta, rep(0, p+1-length(theta)))
public.data <- generate.data(1000, p, theta)
public.data$X <- public.data$X[,1:4]

gamma.grid <- 20*(1:31)-19
wub <- 25

# comparison.data <- get.comparison.data(epsilons, sample.sizes, public.data, wub,
#                                        gamma.grid)
# write.csv(comparison.data, "RegConstData/different_theta.csv", row.names=FALSE)
comparison.data <- read.csv("RegConstData/different_theta.csv")
plot.comparison.data(epsilons, sample.sizes, comparison.data,
                     subfolder="DifferentTheta/")

######## PUBLIC DATA WITH LARGER THETA #######
set.seed(42)
p <- 10
theta <- c(1, -0.5, 0.5, 1, 1.5, -2.5, -2)
theta <- c(theta, rep(0, p+1-length(theta)))
public.data <- generate.data(1000, p, theta)
public.data$X <- public.data$X[,1:6]

gamma.grid <- 20*(1:31)-19
wub <- 36

# comparison.data <- get.comparison.data(epsilons, sample.sizes, public.data, wub,
#                                        gamma.grid)
# write.csv(comparison.data, "RegConstData/larger_theta.csv", row.names=FALSE)
comparison.data <- read.csv("RegConstData/larger_theta.csv")
plot.comparison.data(epsilons, sample.sizes, comparison.data,
                     subfolder="LargerTheta/")

######## PUBLIC DATA WITH X FROM TRUNCNORM DISTRIBUTION #######
set.seed(42)
p <- 10
theta <- c(1, 1, 1, -1.8, -2.2)
theta <- c(theta, rep(0, p+1-length(theta)))
X.dist <- 'truncnorm'
norm.mus <- c(0.25, 0.5, 0.75, 0.5, rep(0.5, p+1-length(theta)))
norm.sds <- 0.3
public.data <- generate.data(1000, p, theta, X.dist=X.dist,
                             norm.mus=norm.mus, norm.sds=norm.sds)
public.data$X <- public.data$X[,1:4]

gamma.grid <- 20*(1:31)-19
wub <- 24

# comparison.data <- get.comparison.data(epsilons, sample.sizes, public.data, wub,
#                                        gamma.grid)
# write.csv(comparison.data, "RegConstData/truncnorm.csv", row.names=FALSE)
comparison.data <- read.csv("RegConstData/truncnorm.csv")
plot.comparison.data(epsilons, sample.sizes, comparison.data,
                     subfolder="TruncNorm/")

# ###### PUBLIC DATA WITH X FROM TRUNCNORM DISTRIBUTION AND DIFFERENT THETA #######
# set.seed(42)
# p <- 10
# theta <- c(1, -1, -1, 1.5, -1.5)
# theta <- c(theta, rep(0, p+1-length(theta)))
# X.dist <- 'truncnorm'
# norm.mus <- c(0.25, 0.5, 0.75, 0.5, rep(0.5, p+1-length(theta)))
# norm.sds <- 0.3
# public.data <- generate.data(1000, p, theta, X.dist=X.dist,
#                              norm.mus=norm.mus, norm.sds=norm.sds)
# public.data$X <- public.data$X[,1:4]
#
# gamma.grid <- 20*(1:31)-19
# wub <- 20
#
# # comparison.data <- get.comparison.data(epsilons, sample.sizes, public.data, wub,
# #                                        gamma.grid)
# # write.csv(comparison.data, "RegConstData/truncnorm_diff_theta.csv", row.names=FALSE)
# comparison.data <- read.csv("RegConstData/truncnorm_diff_theta.csv")
# plot.comparison.data(epsilons, sample.sizes, comparison.data,
#                      subfolder="TruncNormDiffTheta/")

######## PUBLIC DATA WITH p=2 #######
set.seed(42)
p <- 10
theta <- c(1, 1, -1.8)
theta <- c(theta, rep(0, p+1-length(theta)))
public.data <- generate.data(1000, p, theta)
public.data$X <- public.data$X[,1:2]

gamma.grid <- 20*(1:31)-19
wub <- 30

# comparison.data <- get.comparison.data(epsilons, sample.sizes, public.data, wub,
#                                        gamma.grid)
# write.csv(comparison.data, "RegConstData/different_p2.csv", row.names=FALSE)
comparison.data <- read.csv("RegConstData/different_p2.csv")
plot.comparison.data(epsilons, sample.sizes, comparison.data,
                     subfolder="DifferentP2/")

######## PUBLIC DATA WITH p=6 #######
set.seed(42)
p <- 10
theta <- c(1, 1, 1, -1.8, -2.2, 0.9, -1.5)
theta <- c(theta, rep(0, p+1-length(theta)))
public.data <- generate.data(1000, p, theta)
public.data$X <- public.data$X[,1:6]

gamma.grid <- 20*(1:31)-19
wub <- 30

# comparison.data <- get.comparison.data(epsilons, sample.sizes, public.data, wub,
#                                        gamma.grid)
# write.csv(comparison.data, "RegConstData/different_p6.csv", row.names=FALSE)
comparison.data <- read.csv("RegConstData/different_p6.csv")
plot.comparison.data(epsilons, sample.sizes, comparison.data,
                     subfolder="DifferentP6/")

######## PUBLIC DATA WITH p=8 #######
set.seed(42)
p <- 10
theta <- c(1, 1, 1, -1.8, -2.2, 0.9, -1.5, 2, -1)
theta <- c(theta, rep(0, p+1-length(theta)))
public.data <- generate.data(1000, p, theta)
public.data$X <- public.data$X[,1:8]

gamma.grid <- 20*(1:31)-19
wub <- 30

# comparison.data <- get.comparison.data(epsilons, sample.sizes, public.data, wub,
#                                        gamma.grid)
# write.csv(comparison.data, "RegConstData/different_p8.csv", row.names=FALSE)
comparison.data <- read.csv("RegConstData/different_p8.csv")
plot.comparison.data(epsilons, sample.sizes, comparison.data,
                     subfolder="DifferentP8/")

