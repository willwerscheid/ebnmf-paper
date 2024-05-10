# This is my first (rough) analysis of the faces data. If the results
# look compelling, I will put this together into a workflowr analysis.
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(NNLM)
library(flashier)
source("code/faces_functions.R")
load("data/faces.RData")
set.seed(1)

faces_train <- 1 - faces_train
faces_test  <- 1 - faces_test

# It is interesting that the pixel values are very evenly distributed
# between 0 and 1.
hist(faces_train,n = 64)
hist(faces_test,n = 64)

# Plot a sampling of the faces.
n <- nrow(faces_train)
m <- ncol(faces_train)
i <- sample(n,49)
print(plot_faces(faces_train[,i]))

# First learn a "baseline" factor using flashier.
fl_baseline <- flash(faces_train,greedy_Kmax = 1,ebnm_fn = ebnm_normal,
                     var_type = 0,verbose = 0)

# NMF initialized with "baseline" factor.
# Lee & Seung (2001) used K = 49.
set.seed(1)
k <- 49
W0 <- fl_baseline$L_pm
nmf <- nnmf(faces_train,k = k,init = list(W0 = W0),
            method = "scd",max.iter = 100,rel.tol = 1e-8,
            n.threads = 4,verbose = 2)
print(plot_faces(nmf$W))

# flashier with sparse point-exponential prior.
set.seed(1)
nmf0 <- nnmf(faces_train,k = k,init = list(W0 = W0),
             method = "scd",max.iter = 10,rel.tol = 1e-8,
             n.threads = 4,verbose = 2)
nmf0$W <- nmf0$W + 0.0001
nmf0$H <- nmf0$H + 0.0001
pe_prior <- ebnm_point_exponential(x = rep(1,100))
pe_prior$fitted_g$scale <- c(0,1)
pe_prior$fitted_g$pi <- c(0.9999,0.0001)
ebnm_pe_prior <- flash_ebnm(prior_family = "point_exponential",
                            fix_g = TRUE,g_init = pe_prior)
fl_pe <- flash_init(faces_train,var_type = 0,S = 0.01)
fl_pe <- flash_factors_init(fl_pe,
                            list(nmf0$W[,k + 1,drop = FALSE],
                                 t(nmf0$H[k + 1,,drop = FALSE])),
                            ebnm_normal)
fl_pe <- flash_factors_init(fl_pe,
                            list(nmf0$W[,1:k],
                                 t(nmf0$H[1:k,])),
                            ebnm_pe_prior)
fl_pe <- flash_backfit(fl_pe,maxiter = 100,verbose = 3)
L <- ldf(fl_pe,type = "i")$L
print(plot_faces(L))

# flashier with sparse truncated normal prior.
tn_prior <- ebnm_generalized_binary(x = rep(1,100),mode = 0.01,scale = 100)
tn_prior$fitted_g$pi <- c(0.99,0.01)
ebnm_tn_prior <- flash_ebnm(prior_family = "generalized_binary",
                            fix_g = TRUE,g_init = tn_prior)
fl_tn <- flash_init(faces_train,var_type = 0,S = 0.01)
fl_tn <- flash_factors_init(fl_tn,
                            list(nmf0$W[,k + 1,drop = FALSE],
                                 t(nmf0$H[k + 1,,drop = FALSE])),
                            ebnm_normal)
fl_tn <- flash_factors_init(fl_tn,
                            list(nmf0$W[,1:k],
                                 t(nmf0$H[1:k,])),
                            ebnm_tn_prior)
fl_tn <- flash_backfit(fl_tn,maxiter = 100,verbose = 3)
L <- ldf(fl_tn,type = "i")$L
print(plot_faces(L))

stop()

fl_pe <- flash(faces_train,greedy_Kmax = 1,
                    ebnm_fn = ebnm_normal)
nmf0 <- nnmf(faces_train,k = 48,method = "scd",
             max.iter = 4,rel.tol = 1e-8,
             n.threads = 4,verbose = 2)
sparse_prior <- ebnm_point_exponential(x = rep(1,100))
sparse_prior$fitted_g$scale <- c(0,1)
sparse_prior$fitted_g$pi <- c(0.9999,0.0001)
ebnm_sparse_prior <- flash_ebnm(prior_family = "point_exponential",
                                fix_g = TRUE,g_init = sparse_prior)
fit_sparse <- flash_factors_init(fit_sparse,
                                 list(nmf0$W,t(nmf0$H)),
                                 c(ebnm_sparse_prior,ebnm_point_exponential))
fit_sparse <- flash_backfit(fit_sparse,maxiter = 1000,verbose = 2)

W <- normalize.cols(nmf$W)
L <- ldf(fit_sparse,type = "i")$L
L[L < 0] <- 0
L <- normalize.cols(L)
i <- order(colSums(W),decreasing = TRUE)
j <- order(colSums(L),decreasing = TRUE)
print(plot_faces(W[,i]))
print(plot_faces(L[,j]))

x <- quantile(W,seq(0,1,0.01))
y <- quantile(L[,-1],seq(0,1,0.01))
plot(x,y,pch = 20,xlab = "NMF",ylab = "flashier")
abline(a = 0,b = 1,col = "magenta",lty = "dotted")

# TO DO NEXT: Try NMF with Hoyer sparsity penalty.
