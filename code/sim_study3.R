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

ns <- c(rep(400, 4), rep(25, 7))
p <- 2000

sim_data <- function(ns, p, dispersion, n_anchor_words = 3, link = "log1p") {
  pops <- rep(LETTERS[1:length(ns)], times = ns)

  # Loadings (document-topics):
  L <- matrix(0, nrow = sum(ns), ncol = 7)
  L[, 7] <- c(rep(1/3, sum(ns[1:4])), rep(1, ns[5]), rep(0, sum(ns[6:11]))) # root
  L[, 5] <- c(rep(1/3, sum(ns[1:2])), rep(0, sum(ns[3:5])), rep(1, ns[6]), rep(0, sum(ns[7:11]))) # branch 1
  L[, 6] <- c(rep(0, sum(ns[1:2])), rep(1/3, sum(ns[3:4])), rep(0, sum(ns[5:6])), rep(1, ns[7]), rep(0, sum(ns[8:11]))) # branch 2
  L[, 1] <- c(rep(1/3, ns[1]), rep(0, sum(ns[2:7])), rep(1, ns[8]), rep(0, sum(ns[9:11]))) # leaf 1
  L[, 2] <- c(rep(0, ns[1]), rep(1/3, ns[2]), rep(0, sum(ns[3:8])), rep(1, ns[9]), rep(0, sum(ns[10:11]))) # leaf 2
  L[, 3] <- c(rep(0, sum(ns[1:2])), rep(1/3, ns[3]), rep(0, sum(ns[4:9])), rep(1, ns[10]), rep(0, ns[11])) # leaf 3
  L[, 4] <- c(rep(0, sum(ns[1:3])), rep(1/3, ns[4]), rep(0, sum(ns[5:10])), rep(1, ns[11])) # leaf 4

  # vary sizes:
  L <- L * rgamma(sum(ns), shape = 2, rate = 2)

  F <- sim_F(p, 7, gamma_shape = 2/3, gamma_scale = 1, n_anchor_words)
  X <- sim_X(L, F, dispersion = dispersion)
  return_sim_data(X, L, F, pops)
}

ss3_ebnmf_K7_res <- run_ebnmf_sims("Y", 7, verbose = TRUE)
ss3_ebnmf_K10_res <- run_ebnmf_sims("Y", 10, verbose = TRUE)
save(ss3_ebnmf_K7_res, ss3_ebnmf_K10_res, file = "output/ss3_ebnmf.RData")

ss3_rcppML_K7_res <- run_RcppML_sims("Y", 7, verbose = TRUE)
ss3_rcppML_K10_res <- run_RcppML_sims("Y", 10, verbose = TRUE)
save(ss3_rcppML_K7_res, ss3_rcppML_K10_res, file = "ss3_rcppML.RData")

setwd("./matlab/")
Matlab$startServer()
matlab <- Matlab()
open(matlab)
setVerbose(matlab, threshold = 100000)
ss3_hoyer_K7_res <- run_Hoyer_sims(matlab, "Y", 7, verbose = TRUE)
ss3_hoyer_K10_res <- run_Hoyer_sims(matlab, "Y", 10, verbose = TRUE)
close(matlab)
setwd("../")
save(ss3_hoyer_K7_res, ss3_hoyer_K10_res, file = "output/ss3_hoyer.RData")

ss3_nnlm_K7_res <- run_nnlm_sims("Y", 7, verbose = TRUE)
ss3_nnlm_K10_res <- run_nnlm_sims("Y", 10, verbose = TRUE)
save(ss3_nnlm_K7_res, ss3_nnlm_K10_res, file = "output/ss3_nnlm.RData")


### Plots -----

rel_widths <- c(0.6, 0.4)

xlab <- "L1 penalty"

make_trueK_plot(
  ss3_rcppML_K7_res, ss3_ebnmf_K7_res, "RcppML", xlab, ylab, TRUE, c(0.05, 0.25, 0.5)
)
ggsave("output/plots/ss3_rcppML_K7.pdf", width = 10, height = 6)

make_overK_plot(
  ss3_rcppML_K10_res, ss3_ebnmf_K10_res, "RcppML", xlab, ylab, rel_widths, 10, TRUE, c(0.05, 0.25, 0.5)
)
ggsave("output/plots/ss3_rcppML_K10.pdf", width = 12, height = 8)

make_trueK_plot(
  ss3_nnlm_K7_res, ss3_ebnmf_K7_res, "NNLM", xlab, ylab, TRUE, c(0.5, 2, 5)
)
ggsave("output/plots/ss3_nnlm_K7.pdf", width = 10, height = 6)

make_overK_plot(
  ss3_nnlm_K10_res, ss3_ebnmf_K10_res, "NNLM", xlab, ylab, rel_widths, 10, TRUE, c(0.5, 2, 5)
)
ggsave("output/plots/ss3_nnlm_K10.pdf", width = 12, height = 8)

xlab <- "L1 penalty"

make_trueK_plot(
  ss3_hoyer_K7_res, ss3_ebnmf_K7_res, "Hoyer", "Sparseness penalty", ylab
)
ggsave("output/plots/ss3_hoyer_K7.pdf", width = 10, height = 6)

make_overK_plot(
  ss3_hoyer_K10_res, ss3_ebnmf_K10_res, "Hoyer", xlab, ylab, rel_widths, 10, FALSE, c(0, 0.5, 1)
)
ggsave("output/plots/ss3_hoyer_K10.pdf", width = 12, height = 8)
