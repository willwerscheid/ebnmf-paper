# Short script to illustrate the use of SVD and NMF to initialize
# flashier fits.
library(rsvd)
library(NNLM)
library(flashier)
source("code/faces_functions.R")
load("data/faces.RData")
set.seed(1)
faces <- 1 - faces_train
n <- nrow(faces)
m <- ncol(faces)

# EBNMF with random initialization.
k <- 49
fit1 <- flash_init(faces,var_type = 0)
fit1 <- flash_factors_init(fit1,
                           list(matrix(runif(n*k),n,k),
                                matrix(runif(m*k),m,k)),
                           ebnm_point_exponential)
capture.output(fit1 <- flash_backfit(fit1,maxiter = 100,verbose = 0))
print(fit1$elbo)

# EBNMF initialized using NNLM.
set.seed(1)
nmf <- nnmf(faces,k = 49,method = "scd",
            max.iter = 10,rel.tol = 1e-8,
            n.threads = 4,verbose = 0)
fit2 <- flash_init(faces,var_type = 0)
fit2 <- flash_factors_init(fit2,
                           list(nmf$W,t(nmf$H)),
                           ebnm_point_exponential)
fit2 <- flash_backfit(fit2,maxiter = 100,verbose = 0)
print(fit2$elbo)

# EBMF with random initialization.
k <- 24
fit3 <- flash_init(faces,var_type = 0)
fit3 <- flash_factors_init(fit3,
                           list(matrix(rnorm(n*k),n,k),
                                matrix(rnorm(m*k),m,k)),
                           ebnm_point_normal)
capture.output(fit3 <- flash_backfit(fit3,maxiter = 100,verbose = 0))
print(fit3$elbo)

# EBMF initialized with PCA.
res <- rsvd(faces,k = k)
fit4 <- flash_init(faces,var_type = 0)
fit4 <- flash_factors_init(fit4,
                           with(res,
                                list(u %*% diag(sqrt(d)),
                                     v %*% diag(sqrt(d)))),
                           ebnm_point_normal)
fit4 <- flash_backfit(fit4,maxiter = 100,verbose = 0)
print(fit4$elbo)
