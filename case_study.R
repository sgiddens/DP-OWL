####################################################################################
## install the DPpack package
# library(devtools)
# library(roxygen2)
# 
# devtools::document(setwd("~/Library/CloudStorage/OneDrive-St.JudeChildren'sResearchHospital/DP-ITR/real_data_analysis/DPpack-main"))
# devtools::install()
####################################################################################
library(DPpack)

####################################################################################
data_dp_itr <- read.csv("data_dp_itr.txt", header = TRUE)
colnames(data_dp_itr)

X <- data_dp_itr[,c("Age_at_Dx","Age_at_Enrollment","Sex","Race","Diagnosis")]
summary(X$Age_at_Dx)
summary(X$Age_at_Enrollment)
# X$Age_at_Dx <- scale(X$Age_at_Dx, center = TRUE, scale = FALSE)
# X$Age_at_Dx <- (X$Age_at_Dx-min(X$Age_at_Dx))/(max(X$Age_at_Dx)-min(X$Age_at_Dx))
# X$Age_at_Enrollment <- scale(X$Age_at_Enrollment, center = TRUE, scale = FALSE)
# X$Age_at_Enrollment <- (X$Age_at_Enrollment-min(X$Age_at_Enrollment))/(max(X$Age_at_Enrollment)-min(X$Age_at_Enrollment))
X$Sex <- ifelse(X$Sex=="Female", 0, 1)  ##0 (Female): 126; 1 (Male): 120.
X$Race <- ifelse(X$Race=="White", 0, 1)  ##0 (White): 217; 1 (Other): 29.
X$Diagnosis <- ifelse(X$Diagnosis=="Leukemia",0,1)  ##0 (Leukemia): 109; 1 (non-Leukemia): 137.
colnames(X) <- c("Age_at_Dx","Age_at_Enrollment","Sex","Race","Diagnosis")

A <- data_dp_itr$Assigned_Treatment  ##1: placebo (126); 2: Melatonin (120).
A <- ifelse(A==1, -1, 1)  ##-1: placebo (126); 1: Melatonin (120).

B <- data_dp_itr$WASI_Matrix_Reasoning_IQ_Difference_.FU_BL_
hist(B)
summary(B)

B <- B+abs(min(B))+0.001

upper.bounds <- c(20,60,1,1,1)
lower.bounds <- c(0,18,0,0,0)

training.data <- list("X"=X,
                      "A"=A,
                      "B"=B,
                      "upper.bounds"=upper.bounds,
                      "lower.bounds"=lower.bounds)

####################################################################################
kernel <- 'linear'
wub <- 8 # Global/public weights upper bound (sup{benefit/propensity})

####################################################################################
data_par <- read.csv("Regularization constant data.csv", header = TRUE)
unique(data_par$eps)
unique(data_par$gamma[data_par$sample.size==200])

eps <- c(0.1,0.5,1,2,5,Inf)
gamma <- c(30,50,90,100,150)

par_list <- as.data.frame(as.matrix(expand.grid(eps, gamma)))
colnames(par_list) <- c("eps","gamma")

real_data_analysis <- function(i, seed){
  set.seed(seed)
  owl.model <- train.model(training.data = training.data,
                           eps = par_list$eps[i],
                           gamma = par_list$gamma[i],
                           kernel = kernel,
                           weights.upper.bound = wub)
  
  results <- evaluate.model(trained.model = owl.model, 
                            validation.data = training.data,
                            ground.truth = FALSE) 
  
  pred.label <- owl.model$predict(training.data$X, add.bias = TRUE)
  
  return(list("coeff"=owl.model$coeff,
              "evf"=results$val.func,
              "pred.label"=pred.label))
}

####################################################################################
rep <- 200

results_list <- lapply(1:nrow(par_list), function(i) lapply(1:rep, function(j) real_data_analysis(i,seed=1234+j)))

saveRDS(results_list, "results_list.rds") ## save for further results summary
