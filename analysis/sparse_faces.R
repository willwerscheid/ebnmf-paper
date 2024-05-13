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
  nnmf(faces_train,k = k,init = list(W0 = W0),method = "scd",max.iter = 4,
       verbose = 0)
nmf0 <- nnmf(faces_train,k = k + 1,init = list(W = nmf_fixed_baseline$W),
             method = "scd",max.iter = 20,verbose = 0)
nmf <- nnmf(faces_train,k = k,method = "scd",max.iter = 200,rel.tol = 1e-8,
            init = list(W0 = nmf0$W[,k + 1]),n.threads = 4,verbose = 2)

# flashier refinement of NMF solution.
fl_nmf <- flash_init(faces_train,var_type = 0)
fl_nmf <- flash_factors_init(fl_nmf,
                             list(nmf$W,t(nmf$H)),
                             ebnm_point_exponential)
fl_nmf <- flash_factors_fix(fl_nmf,kset = k + 1,which_dim = "loadings")
fl_nmf <- flash_backfit(fl_nmf,maxiter = 100,extrapolate = FALSE,verbose = 2)
fl_nmf <- flash_backfit(fl_nmf,maxiter = 100,extrapolate = TRUE,verbose = 3)
print(plot_faces(fl_nmf$L_pm,title = "NMF",order_by_sparsity = TRUE))

# flashier initialized to NMF, but with much more wiggle room.
set.seed(1)
W0 <- nmf$W[,k + 1]
nmf0 <- nnmf(faces_train,k = k,init = list(W0 = W0),method = "scd",
             max.iter = 10,rel.tol = 1e-8)
fl <- flash_init(faces_train,var_type = 0)
fl <- flash_factors_init(fl,list(nmf0$W,t(nmf0$H)),ebnm_point_exponential)
fl <- flash_factors_fix(fl,kset = k + 1,which_dim = "loadings")
fl <- flash_backfit(fl,maxiter = 100,extrapolate = FALSE,verbose = 2)
fl <- flash_backfit(fl,maxiter = 100,extrapolate = TRUE,verbose = 3)
print(plot_faces(fl$L_pm,title = "flashier",order_by_sparsity = TRUE))

# Compare the ELBOs.
cat(sprintf("NMF:      %0.4f\n",fl_nmf$elbo))
cat(sprintf("flashier: %0.4f\n",fl$elbo))

# Compare sparsity of NMF and flashier/NMF solutions.
zero <- 0.001
pdat <-
  data.frame(nmf      = sort(colMeans(normalize.cols(nmf$W) > zero)),
             flashier = sort(colMeans(normalize.cols(fl_nmf$L_pm) > zero)))
p1 <- ggplot(pdat,aes(x = nmf,y = flashier)) +
  geom_point() +
  geom_abline(intercept = 0,slope = 1,color = "magenta",lty = "dotted") +
  theme_cowplot(font_size = 12)
print(p1)

# Compare sparsity of NMF and flashier solutions.
pdat <-
  data.frame(nmf      = sort(colMeans(normalize.cols(nmf$W) > zero)),
             flashier = sort(colMeans(normalize.cols(fl$L_pm) > zero)))
p2 <- ggplot(pdat,aes(x = nmf,y = flashier)) +
  geom_point() +
  geom_abline(intercept = 0,slope = 1,color = "magenta",lty = "dotted") +
  theme_cowplot(font_size = 12)

# Compare the (Frobenius) NMF objective.
set.seed(1)
fl_nmf_train <- nnmf(faces_train,k = k + 1,
                  init = list(W0 = fl_nmf$L_pm),
                  method = "scd",max.iter = 40,n.threads = 4,verbose = 0)
fl_train <- nnmf(faces_train,k = k + 1,
                 init = list(W0 = fl$L_pm),
                 method = "scd",max.iter = 40,n.threads = 4,verbose = 0)
print(frobenius_norm(faces_train,fl_nmf_train$W,fl_nmf_train$H))
print(frobenius_norm(faces_train,fl_train$W,fl_train$H))
