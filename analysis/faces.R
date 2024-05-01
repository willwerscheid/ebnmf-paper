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

# flashier with fixed prior.
k <- 49
fixed_prior <- ebnm_point_exponential(x = rep(1,100))
fixed_prior$fitted_g$scale <- c(0,1)
fixed_prior$fitted_g$pi <- c(0.1,0.9)
ebnm_fixed_prior <- flash_ebnm(prior_family = "point_exponential",
                               fix_g = TRUE,g_init = fixed_prior)
fit_fixed <- flash_init(faces_train,var_type = 0)
fit_fixed <- flash_factors_init(fit_fixed,
                                list(nmf$W,t(nmf$H)),
                                ebnm_fixed_prior)
fit_fixed <- flash_backfit(fit_fixed,verbose = 2)
print(plot_faces(ldf(fit_fixed,type = "i")$L))

# flashier with adaptive prior.
fit_adapt <- flash_init(faces_train,var_type = 0)
fit_adapt <- flash_factors_init(fit_adapt,
                                list(fit_fixed$L_pm,fit_fixed$F_pm),
                                ebnm_point_exponential)
fit_adapt <- flash_backfit(fit_adapt,verbose = 2)
print(plot_faces(ldf(fit_adaptive,type = "i")$L))
