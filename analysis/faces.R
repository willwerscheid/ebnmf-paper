# This is my first (rough) analysis of the faces data. If the results
# look interesting, I will put this together into a workflowr analysis.
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
i <- sample(n,49)
print(plot_faces(faces_train[,i]))

# TO DO:
#
# Try with k = 98.
#
# Use some of the test data for training.
#

# Lee & Seung (2001) used K = 49.
set.seed(1)
nmf <- nnmf(faces_train,k = 49,method = "scd",
            max.iter = 200,rel.tol = 1e-8,
            n.threads = 4,verbose = 2)
print(plot_faces(nmf$W))

stop()

# flashier with "flat" prior.
k <- 49
set.seed(1)
nmf0 <- nnmf(faces_train,k = 49,method = "scd",
            max.iter = 10,rel.tol = 1e-8,
            n.threads = 4,verbose = 2)
flat_prior <- ebnm_point_exponential(x = rep(1,100))
flat_prior$fitted_g$scale <- c(0,1)
flat_prior$fitted_g$pi <- c(0.9999,0.0001)
ebnm_flat_prior <- flash_ebnm(prior_family = "point_exponential",
                              fix_g = TRUE,g_init = flat_prior)
fit_flat <- flash_init(faces_train,var_type = 0)
fit_flat <- flash_factors_init(fit_flat,
                               list(nmf0$W,t(nmf0$H)),
                               c(ebnm_flat_prior,ebnm_point_exponential))
fit_flat <- flash_backfit(fit_flat,maxiter = 200,verbose = 2)

plot(normalize.cols(nmf$W),normalize.cols(fit_flat$L_pm),pch = 20)

W <- normalize.cols(nmf$W)
L <- normalize.cols(ldf(fit_flat,type = "i")$L)
print(plot_faces(W))
print(plot_faces(L))

# flashier with adaptive prior.
# set.seed(1)
# nmf_init <- nnmf(faces_train,k = 64,method = "scd",max.iter = 10)
fit_adapt <- flash_init(faces_train,var_type = 0)
fit_adapt <- flash_factors_init(fit_adapt,
                                list(nmf$W,t(nmf$H)),
                              # list(fit_fixed$L_pm,fit_fixed$F_pm),
                                ebnm_point_exponential)
fit_adapt <- flash_backfit(fit_adapt,maxiter = 200,verbose = 2)

L  <- ldf(fit_adapt,type = "i")$L
p1 <- plot_faces(nmf$W,nrow = 8,ncol = 8) + ggtitle("NMF")
p2 <- plot_faces(L,nrow = 8,ncol = 8) + ggtitle("flashier")

# Project the training samples.
H_nmf   <- nnlm(nmf$W,faces_train,method = "scd")$coefficients
H_adapt <- nnlm(L,faces_train,method = "scd")$coefficients

# Compute the Frobenius norm on the training data.
cat("nmf:     ",frobenius_norm(faces_train,nmf$W,H_nmf),"\n")
cat("flashier:",frobenius_norm(faces_train,L,H_adapt),"\n")

# Project the test samples.
H_nmf   <- nnlm(nmf$W,faces_test,method = "scd")$coefficients
H_adapt <- nnlm(L,faces_test,method = "scd")$coefficients

# Compute the Frobenius norm on the test data.
cat("nmf:     ",frobenius_norm(faces_test,nmf$W,H_nmf),"\n")
cat("flashier:",frobenius_norm(faces_test,L,H_adapt),"\n")
