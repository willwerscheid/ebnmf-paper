# TO DO: Explain here what this script is for, and how to use it.
library(R.matlab)
library(ebnm)
library(flashier)
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# Load the swimmers data set.
Y <- readMat("../data/swimmer.mat")$Y
Y <- apply(Y,3,as.vector)
Y <- Y - 1

# This reproduces Jason's result:
#
#   fit <- flash(Y,ebnm_fn = ebnm_point_exponential,
#                backfit = TRUE,var_type = 0)
#

fit1 <- flash(Y,greedy_Kmax = 2,
              ebnm_fn = ebnm_point_exponential,
              backfit = TRUE,var_type = 0)
print(plot_images(ldf(fit1)$L))

# Initialize the flash fit.
n <- nrow(Y)
m <- ncol(Y)
k <- 17
gl <- ebnm_point_exponential(x = c(rep(40,100)))
gf <- ebnm_point_exponential(x = c(rep(1,100)))
gl$fitted_g$pi <- c(0.99,0.01)
gf$fitted_g$pi <- c(0.9,0.1)
ebnm_exp_L <- flash_ebnm(prior_family = "point_exponential",
                         fix_g = TRUE,
                         g_init = gl)
ebnm_exp_F <- flash_ebnm(prior_family = "point_exponential",
                         fix_g = TRUE,
                         g_init = gf)
# init <- list(u = matrix(runif(n*k),n,k),
#              v = matrix(runif(m*k),m,k),
#              d = rep(1,256))
fit <- flash_init(Y,var_type = 0)
fit <- flash_factors_init(fit,
                          list(matrix(runif(n*k),n,k),
                               matrix(runif(m*k),m,k)),
                          list(L = ebnm_exp_L,F = ebnm_exp_F))
fit <- flash_backfit(fit)

plot_images <- function (fitted_L) {
  fitted_L <- max(fitted_L) - fitted_L
  colnames(fitted_L) <- paste0("V", 1:ncol(fitted_L))
  tib <- as_tibble(fitted_L)
  tib <- mutate(tib,
    row = rep(1:32,each = 32),
    col = rep(1:32,times = 32)) 
  tib <- pivot_longer(tib,
      cols = -c(row,col),
      names_to = "k",
      values_to = "loading",
      names_prefix = "V",
      names_transform = as.numeric)
  p <- ggplot(tib,aes(x = row,y = col,fill = loading)) +
    geom_tile() +
    scale_y_reverse() +
    scale_fill_gradient(low = "white",high = "black") +
    facet_wrap(~k) +
    guides(fill = "none",x = "none",y = "none") +
    theme_cowplot(font_size = 8)
  return(p)
}

print(plot_images(ldf(fit)$L))
