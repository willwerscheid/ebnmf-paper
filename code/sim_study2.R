library(R.matlab)
library(tibble)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(ebnm)
library(flashier)
library(fastTopics)
source("./code/sim_functions.R")
source("./code/sim_study_functions.R")
options(matlab = "/Applications/MATLAB_R2025a.app/bin/matlab")


### Simulations -----

ns <- c(rep(25, 3), rep(500, 3))
p <- 2000

sim_data <- function(ns, p, dispersion, n_anchor_words = 3) {
  pops <- rep(LETTERS[1:length(ns)], times = ns)

  L <- matrix(0, nrow = sum(ns), ncol = 3)
  pi1 <- seq(0, 1, length.out = ns[4])
  pi2 <- seq(0, 1, length.out = ns[5])
  pi3 <- seq(0, 1, length.out = ns[6])
  L[, 1] <- c(rep(1, ns[1]), rep(0, sum(ns[2:3])), pi1, pi2, rep(0, ns[6]))
  L[, 2] <- c(rep(0, ns[1]), rep(1, ns[2]), rep(0, ns[3]), 1 - pi1, rep(0, ns[5]), pi3)
  L[, 3] <- c(rep(0, sum(ns[1:2])), rep(1, ns[3]), rep(0, ns[4]), 1 - pi2, 1 - pi3)

  # vary sizes:
  L <- L * rgamma(sum(ns), shape = 2, rate = 2)

  F <- sim_F(p, 3, gamma_shape = 2/3, gamma_scale = 1, n_anchor_words)
  X <- sim_X(L, F, dispersion = dispersion)
  return_sim_data(X, L, F, pops)
}

ss2_ebnmf_K3_res <- run_ebnmf_sims("Y", 3, verbose = TRUE)
ss2_ebnmf_K6_res <- run_ebnmf_sims("Y", 6, verbose = TRUE)
save(ss2_ebnmf_K3_res, ss2_ebnmf_K6_res, file = "output/ss2_ebnmf.RData")

ss2_rcppML_K3_res <- run_RcppML_sims("Y", 3, verbose = TRUE)
ss2_rcppML_K6_res <- run_RcppML_sims("Y", 6, verbose = TRUE)
save(ss2_rcppML_K3_res, ss2_rcppML_K6_res, file = "output/ss2_rcppML.RData")

setwd("./matlab/")
Matlab$startServer()
matlab <- Matlab()
open(matlab)
setVerbose(matlab, threshold = 100000)
ss2_hoyer_K3_res <- run_Hoyer_sims(matlab, "Y", 3, verbose = TRUE)
ss2_hoyer_K6_res <- run_Hoyer_sims(matlab, "Y", 6, verbose = TRUE)
close(matlab)
setwd("../")
save(ss2_hoyer_K3_res, ss2_hoyer_K6_res, file = "output/ss2_hoyer.RData")

ss2_nnlm_K3_res <- run_nnlm_sims("Y", 3, verbose = TRUE)
ss2_nnlm_K6_res <- run_nnlm_sims("Y", 6, verbose = TRUE)
save(ss2_nnlm_K3_res, ss2_nnlm_K6_res, file = "output/ss2_nnlm.RData")


### Plots -----

rel_widths <- c(0.5, 0.5)

xlab <- "L1 penalty"

make_trueK_plot(
  ss2_rcppML_K3_res, ss2_ebnmf_K3_res, "RcppML", xlab, ylab, TRUE, c(0.025, 0.1, 0.25, 0.5)
)
ggsave("output/plots/ss2_rcppML_K3.pdf", width = 10, height = 6)

make_overK_plot(
  ss2_rcppML_K6_res, ss2_ebnmf_K6_res, "RcppML", xlab, ylab, rel_widths, 6, TRUE, c(0.05, 0.25, 0.5)
)
ggsave("output/plots/ss2_rcppML_K6.pdf", width = 12, height = 8)

make_trueK_plot(
  ss2_nnlm_K3_res, ss2_ebnmf_K3_res, "NNLM", xlab, ylab, TRUE, c(0.5, 2, 5)
)
ggsave("output/plots/ss2_nnlm_K3.pdf", width = 10, height = 6)

make_overK_plot(
  ss2_nnlm_K6_res, ss2_ebnmf_K6_res, "NNLM", xlab, ylab, rel_widths, 6, TRUE, c(0.5, 2, 5)
)
ggsave("output/plots/ss2_nnlm_K6.pdf", width = 12, height = 8)

xlab <- "L1 penalty"

make_trueK_plot(
  ss2_hoyer_K3_res, ss2_ebnmf_K3_res, "Hoyer", "Sparseness penalty", ylab
)
ggsave("output/plots/ss2_hoyer_K3.pdf", width = 10, height = 6)

make_overK_plot(
  ss2_hoyer_K6_res, ss2_ebnmf_K6_res, "Hoyer", xlab, ylab, rel_widths, 6
)
ggsave("output/plots/ss2_hoyer_K6.pdf", width = 12, height = 8)
