library(RcppML)
library(NNLM)
source("code/faces_functions.R")
load("data/faces.RData")
set.seed(1)
faces <- 1 - faces_train
set.seed(1)

# RcppML.
setRcppMLthreads(4)
t0 <- proc.time()
fit1 <- RcppML::nmf(faces,k = 49,maxit = 150,tol = 1e-8,diag = FALSE,
                    verbose = TRUE)
t1 <- proc.time()
print(t1 - t0)

# NNLM.
t0 <- proc.time()
fit2 <- NNLM::nnmf(faces,k = 49,method = "scd",
                   max.iter = 100,rel.tol = 1e-8,
                   n.threads = 4,verbose = 2)
t1 <- proc.time()
print(t1 - t0)
print(frobenius_norm(faces,fit1$w,fit1$h))
print(frobenius_norm(faces,fit2$W,fit2$H))


