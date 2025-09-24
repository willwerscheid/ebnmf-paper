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

ns <- c(1250, 250, 100, 10)
p <- 2000

sim_data <- function(ns, p, dispersion, n_anchor_words = 3) {
  pops <- rep(LETTERS[1:length(ns)], times = ns)

  L <- matrix(0, nrow = sum(ns), ncol = 4)
  L[, 1] <- c(rep(1, ns[1]), rep(0, sum(ns[2:4])))
  L[, 2] <- c(rep(0, ns[1]), rep(1, ns[2]), rep(0, sum(ns[3:4])))
  L[, 3] <- c(rep(0, sum(ns[1:2])), rep(1, ns[3]), rep(0, ns[4]))
  L[, 4] <- c(rep(0, sum(ns[1:3])), rep(1, ns[4]))

  # vary sizes:
  L <- L * rgamma(sum(ns), shape = 2, rate = 2)

  F <- sim_F(p, 4, gamma_shape = 2/3, gamma_scale = 1, n_anchor_words)
  X <- sim_X(L, F, dispersion = dispersion)
  return_sim_data(X, L, F, pops)
}

ss1_ebnmf_K4_res <- run_ebnmf_sims("Y", 4, verbose = TRUE)
ss1_ebnmf_K6_res <- run_ebnmf_sims("Y", 6, verbose = TRUE)
save(ss1_ebnmf_K4_res, ss1_ebnmf_K6_res, file = "output/ss1_ebnmf.RData")

ss1_rcppML_K4_res <- run_RcppML_sims("Y", 4, verbose = TRUE)
ss1_rcppML_K6_res <- run_RcppML_sims("Y", 6, verbose = TRUE)
save(ss1_rcppML_K4_res, ss1_rcppML_K6_res, file = "ss1_rcppML.RData")

setwd("./matlab/")
Matlab$startServer()
matlab <- Matlab()
open(matlab)
setVerbose(matlab, threshold = 100000)
ss1_hoyer_K4_res <- run_Hoyer_sims(matlab, "Y", 4, verbose = TRUE)
ss1_hoyer_K6_res <- run_Hoyer_sims(matlab, "Y", 6, verbose = TRUE)
close(matlab)
setwd("../")
save(ss1_hoyer_K4_res, ss1_hoyer_K6_res, file = "output/ss1_hoyer.RData")

ss1_nnlm_K4_res <- run_nnlm_sims("Y", 4, verbose = TRUE)
ss1_nnlm_K6_res <- run_nnlm_sims("Y", 6, verbose = TRUE)
save(ss1_nnlm_K4_res, ss1_nnlm_K6_res, file = "output/ss1_nnlm.RData")


### Plots -----

rel_widths <- c(0.6, 0.4)

xlab <- "L1 penalty"

make_trueK_plot(
  ss1_rcppML_K4_res, ss1_ebnmf_K4_res, "RcppML", xlab, ylab, TRUE, c(0.025, 0.1, 0.25, 0.5)
)
ggsave("output/plots/ss1_rcppML_K4.pdf", width = 10, height = 6)

make_overK_plot(
  ss1_rcppML_K6_res, ss1_ebnmf_K6_res, "RcppML", xlab, ylab, rel_widths, 6, TRUE, c(0.05, 0.25, 0.5)
)
ggsave("output/plots/ss1_rcppML_K6.pdf", width = 12, height = 8)

make_trueK_plot(
  ss1_nnlm_K4_res, ss1_ebnmf_K4_res, "NNLM", xlab, ylab, TRUE, c(0.5, 2, 5)
)
ggsave("output/plots/ss1_nnlm_K4.pdf", width = 10, height = 6)

make_overK_plot(
  ss1_nnlm_K6_res, ss1_ebnmf_K6_res, "NNLM", xlab, ylab, rel_widths, 6, TRUE, c(0.5, 2, 5)
)
ggsave("output/plots/ss1_nnlm_K6.pdf", width = 12, height = 8)

xlab <- "L1 penalty"

make_trueK_plot(
  ss1_hoyer_K4_res, ss1_ebnmf_K4_res, "Hoyer", "Sparseness penalty", ylab
)
ggsave("output/plots/ss1_hoyer_K4.pdf", width = 10, height = 6)

make_overK_plot(
  ss1_hoyer_K6_res, ss1_ebnmf_K6_res, "Hoyer", xlab, ylab, rel_widths, 6
)
ggsave("output/plots/ss1_hoyer_K6.pdf", width = 12, height = 8)
