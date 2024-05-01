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
faces_test <- 1 - faces_test

# It is interesting that the pixel values are very evenly distributed
# between 0 and 1.
hist(faces_train,n = 64)
hist(faces_test,n = 64)

# Plot a sampling of the faces.
n <- nrow(faces_train)
i <- sample(n,49)
print(plot_faces(faces_train[,i]))

# Lee & Seung (2001) used K = 49.
nmf <- nnmf(faces_train,k = 49,method = "scd",
            max.iter = 200,rel.tol = 1e-8,
            n.threads = 4,verbose = 2)
print(plot_faces(nmf$W))

# NMF with k = 25.
set.seed(1)
nmf <- nnmf(faces_train,k = 25,method = "scd",
            max.iter = 200,rel.tol = 1e-8,
            n.threads = 4,verbose = 2)

# flashier with fixed prior.
k <- 25
fixed_prior <- ebnm_point_exponential(x = rep(1,100))
fixed_prior$fitted_g$scale <- c(0,1)
fixed_prior$fitted_g$pi <- c(0.1,0.9)
ebnm_fixed_prior <- flash_ebnm(prior_family = "point_exponential",
                               fix_g = TRUE,g_init = fixed_prior)
fit_fixed <- flash_init(faces_train,var_type = 0)
fit_fixed <- flash_factors_init(fit_fixed,
                                list(nmf$W,t(nmf$H)),
                                ebnm_fixed_prior)
fit_fixed <- flash_backfit(fit_fixed,maxiter = 100,verbose = 2)

# flashier with adaptive prior.
fit_adapt <- flash_init(faces_train,var_type = 0)
fit_adapt <- flash_factors_init(fit_adapt,
                                list(fit_fixed$L_pm,fit_fixed$F_pm),
                                ebnm_point_exponential)
fit_adapt <- flash_backfit(fit_adapt,maxiter = 1000,verbose = 2)

p1 <- plot_faces(nmf$W) + ggtitle("NMF")
p2 <- plot_faces(ldf(fit_fixed,type = "i")$L) + ggtitle("flashier_fixed")
p3 <- plot_faces(ldf(fit_adapt,type = "i")$L) + ggtitle("flashier_adapt")
plot_grid(p1,p3,nrow = 2,ncol = 1)

