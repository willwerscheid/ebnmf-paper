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
# Lee & Seung (2001) used k = 49.
set.seed(1)
k <- 49
W0 <- fl_baseline$L_pm
nmf_fixed_baseline <-
  nnmf(faces_train,k = k,init = list(W0 = W0),method = "scd",max.iter = 4)
nmf0 <- nnmf(faces_train,k = k + 1,init = list(W = nmf_fixed_baseline$W),
             method = "scd",max.iter = 20)
nmf <- nnmf(faces_train,k = k,method = "scd",max.iter = 200,rel.tol = 1e-8,
            init = list(W0 = nmf0$W[,k + 1]),
            n.threads = 4,verbose = 2)
print(plot_faces(nmf$W,title = "NMF",order_by_sparsity = TRUE))

# flashier with point-exponential prior.
set.seed(1)
W0 <- nmf$W[,k + 1]
nmf0 <- nnmf(faces_train,k = k,init = list(W0 = W0),method = "scd",
             max.iter = 10,rel.tol = 1e-8)
fl <- flash_init(faces_train,var_type = 0)
fl <- flash_factors_init(fl,
                         list(nmf0$W,t(nmf0$H)),
                         ebnm_point_exponential)
fl <- flash_factors_fix(fl,kset = k + 1,which_dim = "loadings")
fl <- flash_backfit(fl,maxiter = 100,extrapolate = FALSE,verbose = 2)
fl <- flash_backfit(fl,maxiter = 100,extrapolate = TRUE,verbose = 3)
print(plot_faces(fl$L_pm,title = "flashier",order_by_sparsity = TRUE))

# Compare sparsity of NMF and flashier solutions.
plot(sort(colMeans(normalize.cols(nmf$W) > 0.001)),
     sort(colMeans(normalize.cols(fl$L_pm) > 0.001)),
     pch = 20)
abline(a = 0,b = 1,col = "magenta",lty = "dotted")

# TO DO NEXT:
#
#   + Try NMF with Hoyer sparsity penalty.
#  
#   + Assess fit in the test set.
#
