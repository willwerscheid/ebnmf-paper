library(tidyverse)
library(flashier)
source("./code/sim_functions.R")
source("./code/swimmer_functions.R")

matlab_bin <- "/Applications/MATLAB_R2025a.app/bin/matlab"

swimmer <- R.matlab::readMat("./data/swimmer.mat")
dat <- apply(swimmer$Y, 3, as.vector)
dat <- dat - 1

add_noise <- function(dat, noise_lvl) {
  shape1 <- ifelse(dat == 0, noise_lvl, 1)
  shape2 <- ifelse(dat == 0, 1, noise_lvl)
  noisy_dat <- matrix(rbeta(length(dat), shape1, shape2), nrow = nrow(dat))
  return(noisy_dat)
}

# Sample of noisy swimmer images.
set.seed(1)
noisy_dat <- add_noise(dat, 0.1)
plot_images(noisy_dat[, sample(ncol(noisy_dat), 6)], nrow = 1, ncol = 6) +
  theme(strip.text = element_blank())
ggsave("./output/plots/swimmer_sample.pdf", width = 7, height = 1.5)

# "True" factors.
parts <- extract_true_swimmer_factors(dat)
parts <- parts[, c(1, 13, 12,  4, 2, 10, 5, 11,  6,
                      16,  8, 17, 9, 14, 3,  7, 15)]
p_true <- plot_images(parts, nrow = 1, ncol = 17) +
  ggtitle("Ideal parts-based representation ")
ggsave("./output/plots/swimmer_true_factors.pdf", p_true, width = 8, height = 0.75)

# EBNMF results.
fl <- flash(noisy_dat, ebnm_fn = ebnm_point_exponential, greedy_Kmax = 20, backfit = TRUE)
p_ebnmf <- plot_images(fl$L_pm[, align_cols(fl$L_pm, parts)[1:17]],
                       nrow = 1, ncol = 17) +
  ggtitle("EBNMF")
ggsave("./output/plots/swimmer_flash.pdf", p_ebnmf, width = 8, height = 0.9)

# Sparse NMF results (DeBruine).
plotlist <- list()
for (pen in c(0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)) {
  if (pen > 1) {
    penF <- 0.01
  } else {
    penF <- 0
  }
  spnmf <- run_RcppML_sparse_nmf(noisy_dat, k = 20, L1pen = pen, L1pen_F = penF, seeds = 1:10)
  W <- spnmf$fit$W
  p <- plot_images(W[, align_cols(W, parts)[1:17]], nrow = 1, ncol = 17)
  if (pen > 0) {
    p <- p + ggtitle(paste("L1 penalty =", pen))
  } else {
    p <- p + ggtitle("No penalty (vanilla NMF)")
  }
  plotlist <- c(plotlist, list(p))
}
plot_grid(plotlist = plotlist, ncol = 1)
ggsave("./output/plots/swimmer_debruine.pdf", width = 8, height = 8)

setwd("./matlab/")
Matlab$startServer(matlab = matlab_bin)
matlab <- Matlab()
open(matlab)
setVerbose(matlab, threshold = 1000000)

# Sparse NMF results (Hoyer).
plotlist <- list()
for (sprs in c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)) {
  if (sprs > 1) {
    sprsF <- 0.01
  } else {
    sprsF <- 0
  }
  spnmf <- run_Hoyer_sparse_nmf(matlab, noisy_dat, k = 20, sprs = sprs, sprs_F = sprsF, seed = 1:10)
  W <- spnmf$fit$W
  p <- plot_images(W[, align_cols(W, parts)[1:17]], nrow = 1, ncol = 17)
  if (sprs > 0) {
    p <- p + ggtitle(paste("Minimum average sparsity =", sprs))
  } else {
    p <- p + ggtitle("No constraint (vanilla NMF)")
  }
  plotlist <- c(plotlist, list(p))
}
plot_grid(plotlist = plotlist, ncol = 1)
ggsave("../output/plots/swimmer_hoyer.pdf", width = 8, height = 6)

# Sparse NMF results (Kim and Park).
plotlist <- list()
for (pen in c(0, 0.1, 0.3, 1, 3, 10, 30, 100)) {
  if (pen > 0) {
    L2pen <- 1
  } else {
    L2pen <- 0
  }
  spnmf <- run_KimPark_sparse_nmf(matlab, noisy_dat, k = 20, L1pen = pen, L2pen = L2pen, seed = 1:10)
  W <- spnmf$fit$W
  W <- t(t(W) / apply(W, 2, max))
  if (ncol(W) < 17) {
    W <- cbind(
      W, matrix(1e-6 * rexp(nrow(W) * (17 - ncol(W))), nrow = nrow(W), ncol = 17 - ncol(W))
    )
  }
  p <- plot_images(W[, align_cols(W, parts)[1:17]], nrow = 1, ncol = 17, scale_W = FALSE)
  if (pen > 0) {
    p <- p + ggtitle(paste("L1 penalty =", pen))
  } else {
    p <- p + ggtitle("No penalty (vanilla NMF)")
  }
  plotlist <- c(plotlist, list(p))
}
p_kimpark <- plot_grid(plotlist = plotlist, ncol = 1)
ggsave("../output/plots/swimmer_kimpark.pdf", p_kimpark, width = 8, height = 7)

close(matlab)
setwd("../")

plot_grid(plotlist = list(p_true, p_kimpark, p_ebnmf), ncol = 1,
          rel_heights = c(1, 8, 1),
          label_size = 12, labels = "auto")
ggsave("./output/plots/swimmer_combined.pdf", width = 8, height = 9)
