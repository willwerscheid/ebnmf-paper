library(RcppML)
library(NNLM)
frobenius_norm <- function (X, W, H)
  sum((X - W %*% H)^2)
load("data/faces.RData")
faces_train <- 1 - faces_train
set.seed(1)

# Run RcppML method.
set.seed(1)
t0 <- proc.time()
fit1 <- nmf(faces_train,k = 49,maxit = 100,tol = 1e-8,diag = FALSE)
t1 <- proc.time()
print(t1 - t0)

# Run NNLM method.
t0 <- proc.time()
fit2 <- nnmf(faces_train,k = 49,method = "scd",
             max.iter = 100,rel.tol = 1e-8,
             n.threads = 4,verbose = 2)
t1 <- proc.time()
print(t1 - t0)

# Compute the Frobenius norm on the training data.
cat("RcppML:",frobenius_norm(faces_train,fit1$w,fit1$h),"\n")
cat("NNLM:  ",frobenius_norm(faces_train,fit2$W,fit2$H),"\n")
# RcppML: 6704.852
# NNLM:   6596.843
