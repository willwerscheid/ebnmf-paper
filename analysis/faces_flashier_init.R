library(NNLM)
library(flashier)
source("code/faces_functions.R")
load("data/faces.RData")
set.seed(1)
faces <- 1 - faces_train
n <- nrow(faces)
m <- ncol(faces)

# Try with random initialization.
k <- 49
fit1 <- flash_init(faces,var_type = 0)
fit1 <- flash_factors_init(fit1,
                           list(matrix(runif(n*k),n,k),
                                matrix(runif(m*k),m,k)),
                           ebnm_point_exponential)
fit1 <- flash_backfit(fit1,maxiter = 200,verbose = 2)
print(fit1$elbo)
