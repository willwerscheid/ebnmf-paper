library(R.matlab)
library(tibble)
library(dplyr)
library(tidyr)
library(forcats)
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

sim_data <- function(ns, p, dispersion, n_anchor_words = 3) {
  pops <- rep(LETTERS[1:length(ns)], times = ns)

  L <- matrix(0, nrow = sum(ns), ncol = length(ns))
  L[pops == "A", 1] <- 1
  L[pops == "B", 2] <- 1
  L[pops == "C", 3] <- 1

  # vary sizes:
  L <- L * rgamma(sum(ns), shape = 5, scale = 1/5)

  component_scales <- c(1, 1, 1) / 4

  F <- sim_F(p, ncol(L), gamma_shape = 0.5, gamma_scale = 2, n_anchor_words)
  F <- t(t(F) * component_scales)
  X <- sim_X(L, F, dispersion = dispersion)
  return_sim_data(X, L, F, pops)
}

run_sims <- function(matlab, which_dat, Kmax, nreps = 5, verbose = FALSE) {
  all_res <- tibble()
  next_seed <- 0
  for (varied_n in c(2, 3, 5, 10, 15, 20, 25, 30)) {
    for (reps in 1:nreps) {
      if (verbose) cat("RARE N:", varied_n, " SEED:", next_seed + 1, "\n")

      ns <- c(1000, 250, varied_n)
      p <- 2000
      disp <- Inf

      next_seed <- next_seed + 1
      set.seed(next_seed)
      sim_dat <- sim_data(ns, p, disp)

      dat <- sim_dat[[which_dat]]

      if (verbose) cat(" NMF...\n")
      nmf_res <-  run_RcppML_sparse_nmf(dat, k = Kmax, L1pen = 0, seeds = 1:10)
      all_res <- all_res |>
        bind_rows(next_tib(next_seed, disp, varied_n, "NMF", "NMF", Kmax, nmf_res, sim_dat))

      if (verbose) cat(" EBNMF, NMF init...\n")
      ebnmf_res <- run_ebnmf_from_nmf(dat, nmf_res$fit, var_type = 0)
      all_res <- all_res |>
        bind_rows(next_tib(next_seed, disp, varied_n, "EBNMF", "NMF init", Kmax, ebnmf_res, sim_dat))

      if (verbose) cat(" EBNMF, greedy/backfit...\n")
      ebnmf_res <- run_greedy_backfit(dat, Kmax = 6, var_type = 0)
      all_res <- all_res |>
        bind_rows(next_tib(next_seed, disp, varied_n, "EBNMF", "greedy/backfit", Kmax, ebnmf_res, sim_dat))

      if (verbose) cat(" EBNMF, alternating add-many/bf...\n")
      ebnmf_res <- run_greedy_backfit(dat, Kmax = 6, var_type = 0)
      all_res <- all_res |>
        bind_rows(next_tib(next_seed, disp, varied_n, "EBNMF", "alt add-many/bf", Kmax, ebnmf_res, sim_dat))

      if (verbose) cat(" EBNMF, alternating add-one/bf...\n")
      ebnmf_res <- run_alternating(dat, Kmax = 6, var_type = 0)
      all_res <- all_res |>
        bind_rows(next_tib(next_seed, disp, varied_n, "EBNMF", "alt add-one/bf", Kmax, ebnmf_res, sim_dat))

      L1pens2 <- c(1, 3, 10, 30, 100, 300, 1000, 3000, 10000)
      for (L1pen in L1pens2) {
        if (verbose) cat(" Kim Park NMF, L1 =", L1pen, "...\n")
        spnmf_res <- run_KimPark_sparse_nmf(matlab, dat, k = Kmax, L1pen = L1pen, L2pen = 1, seeds = 1:3)
        all_res <- all_res |>
          bind_rows(next_tib(next_seed, disp, varied_n, "KimPark", as.character(L1pen), Kmax, spnmf_res, sim_dat))
      }

      L1pens <- c(0.01, 0.05, seq(0.1, 0.6, by = 0.1))
      for (L1pen in L1pens) {
        if (verbose) cat(" RcppML NMF, L1 =", L1pen, "...\n")
        spnmf_res <- run_RcppML_sparse_nmf(dat, k = Kmax, L1pen = L1pen, seeds = 1:10)
        all_res <- all_res |>
          bind_rows(next_tib(next_seed, disp, varied_n, "RcppML", as.character(L1pen), Kmax, spnmf_res, sim_dat))
      }

      for (L1pen in L1pens) {
        if (verbose) cat(" Hoyer NMF, L1 =", L1pen, "...\n")
        spnmf_res <- run_Hoyer_sparse_nmf(matlab, dat, k = Kmax, pen = L1pen, seeds = 1:3)
        all_res <- all_res |>
          bind_rows(next_tib(next_seed, disp, varied_n, "Hoyer", as.character(L1pen), Kmax, spnmf_res, sim_dat))
      }
    }
  }

  all_res <- all_res |>
    mutate(method = factor(method, levels = unique(all_res$method))) |>
    mutate(submethod = factor(submethod, levels = unique(all_res$submethod)))

  return(all_res)
}

setwd("./matlab/")
Matlab$startServer(matlab = "/Applications/MATLAB_R2025a.app/bin/matlab")
matlab <- Matlab()
open(matlab)
setVerbose(matlab, threshold = 1000000)

sims_K3 <- run_sims(matlab, "Ylibnorm", Kmax = 3, verbose = TRUE)
sims_K4 <- run_sims(matlab, "Ylibnorm", Kmax = 4, verbose = TRUE)
sims_K6 <- run_sims(matlab, "Ylibnorm", Kmax = 6, verbose = TRUE)

close(matlab)
setwd("../")

saveRDS(sims_K3, "./output/ss1_K3.rds")
saveRDS(sims_K4, "./output/ss1_K4.rds")
saveRDS(sims_K6, "./output/ss1_K6.rds")


### Plots -----

plot_df <- sims_K3 |> mutate(K = "K = 3") |>
  bind_rows(sims_K4 |> mutate(K = "K = 4")) |>
  bind_rows(sims_K6 |> mutate(K = "K = 6")) |>
  group_by(K, method, submethod, varied_n, seed) |>
  summarize(metric_val = max(metric_val * (metric_type == "LLcosine3"))) |>
  mutate(metric_val = coalesce(metric_val, 0)) |>
  mutate(varied_n = factor(varied_n)) |>
  mutate(method = fct_rev(method), submethod = fct_rev(submethod)) |>
  group_by(K, method, submethod, varied_n) |>
  summarize(metric_val = mean(metric_val)) |>
  ungroup()

all_plots <- list()
methods <- c("NMF", "EBNMF", "KimPark", "RcppML", "Hoyer")
margins <- c(46, 1, 40, 48, 48)
for (i in 1:length(methods)) {
  next_df <- plot_df |>
    filter(method == methods[[i]]) |>
    mutate(submethod = fct_drop(submethod))
  all_plots[[i]] <- ggplot(next_df, aes(x = varied_n, y = submethod, fill = metric_val)) +
    geom_tile() +
    scale_fill_gradient(low = "red4", high = "seagreen4") +
    theme_minimal() +
    labs(x = "", y = methods[[i]], fill = "") +
    theme(plot.margin = margin(l = margins[i])) +
    facet_wrap(~K, nrow = 1)
  if (methods[[i]] == "NMF") {
    all_plots[[i]] <- all_plots[[i]] +
      labs(y = "")
  } else {
    all_plots[[i]] <- all_plots[[i]] +
      theme(
        strip.background = element_blank(),
        strip.text.x = element_blank()
      )
  }
  if (i < length(methods)) {
    all_plots[[i]] <- all_plots[[i]] +
      theme(axis.text.x = element_blank())
  }
}
legend <- get_legend(all_plots[[1]])
for (i in 1:length(methods)) {
  all_plots[[i]] <- all_plots[[i]] +
    theme(legend.position = "none")
}
p1 <- plot_grid(plotlist = all_plots, ncol = 1, rel_heights = 3 + c(3, 4, 9, 8, 9))
p2 <- plot_grid(p1, legend, nrow = 1, rel_widths = c(9, 1))
p3 <- add_sub(p2, label = "Size of rare population",
              hjust = 0.4, vjust = -0.8, size = 10)
plot(p3)
save_plot("./output/plots/ss1_LL3_allK.pdf", p3, base_width = 7, base_height = 5)
