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
options(matlab = "/Applications/MATLAB_R2025a.app/bin/matlab")

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

run_RcppML_sims <- function(which_dat, Kmax, niter = 5, verbose = FALSE) {
  all_res <- tibble()
  next_seed <- 0
  for (iter in 1:niter) {
    next_seed <- next_seed + 1
    set.seed(next_seed)

    for (disp in c(1, Inf)) {
      if (verbose) cat("ITER: ", iter, "DISPERSION: ", disp, "\n")

      sim_dat <- sim_data(ns, p, disp)
      dat <- sim_dat[[which_dat]]

      L1pens <- seq(0, 0.9, by = 0.05)^2
      for (L1pen in L1pens) {
        if (verbose) cat(" L1 = ", L1pen, "...\n")
        spnmf_res <- run_RcppML_sparse_nmf(dat, k = Kmax, L1pen = L1pen, seeds = 1:10)
        all_res <- all_res |>
          bind_rows(next_tib(next_seed, disp, min(ns), L1pen, Kmax, spnmf_res, sim_dat))
      }
    }
  }

  return(all_res)
}

run_nnlm_sims <- function(which_dat, Kmax, niter = 3, verbose = FALSE) {
  all_res <- tibble()
  next_seed <- 0
  for (iter in 1:niter) {
    next_seed <- next_seed + 1
    set.seed(next_seed)

    for (disp in c(1, Inf)) {
      if (verbose) cat("ITER: ", iter, "DISPERSION: ", disp, "\n")

      sim_dat <- sim_data(ns, p, disp)
      dat <- sim_dat[[which_dat]]

      L1pens <- seq(0, sqrt(10), length.out = 10)^2
      for (L1pen in L1pens) {
        if (verbose) cat(" L1 = ", L1pen, "...\n")
        spnmf_res <- run_sparse_nmf(dat, k = Kmax, L1pen = L1pen, seeds = 1)
        all_res <- all_res |>
          bind_rows(next_tib(next_seed, disp, min(ns), L1pen, Kmax, spnmf_res, sim_dat))
      }
    }
  }

  return(all_res)
}

run_Hoyer_sims <- function(matlab, which_dat, Kmax, niter = 5, verbose = FALSE) {
  all_res <- tibble()
  next_seed <- 0
  for (iter in 1:niter) {
    next_seed <- next_seed + 1
    set.seed(next_seed)

    for (disp in c(1, Inf)) {
      if (verbose) cat("RARE N: ", 10, "DISPERSION: ", disp, "\n")

      sim_dat <- sim_data(ns, p, disp)
      dat <- sim_dat[[which_dat]]

      pens <- seq(0, 0.9, by = 0.05)
      for (pen in pens) {
        if (verbose) cat(" penalty = ", pen, "...\n")

        setVariable(matlab, dat = dat)
        evaluate(matlab, paste0("k=", Kmax,";"))
        evaluate(matlab, paste0("options.sW=", pen,";"))
        evaluate(matlab, paste0("options.sH=0.1;"))
        evaluate(matlab, paste0("options.maxiter=200;"))
        evaluate(matlab, paste0("options.delta=1e-3;"))

        ntrials <- 5
        all_t <- 0
        all_mse <- numeric(ntrials)
        best_mse <- Inf

        for (seed in 1:ntrials) {
          if (verbose) cat("  trial = ", seed, "...\n")

          evaluate(matlab, paste0("rng(", seed, ");"))
          evaluate(matlab, "[W,H,e,t] = sparseNMF(dat,k,options);")

          hoyer_res <- getVariable(matlab, c("W", "H", "e", "t"))
          all_t <- all_t + max(hoyer_res$t)
          all_mse[seed] <- min(hoyer_res$e)
          if (min(hoyer_res$e) < best_mse) {
            best_mse <- min(hoyer_res$e)
            best_fit <- list(W = hoyer_res$W, H = hoyer_res$H)
          }
        }
        spnmf_res <- list(
          t = all_t,
          fit = best_fit,
          all_mse = all_mse
        )
        all_res <- all_res |>
          bind_rows(next_tib(next_seed, disp, min(ns), pen, Kmax, spnmf_res, sim_dat))
      }
    }
  }

  return(all_res)
}

run_ebnmf_sims <- function(which_dat, Kmax, niter = 5, verbose = FALSE) {
  all_res <- tibble()
  next_seed <- 0
  for (iter in 1:niter) {
    next_seed <- next_seed + 1
    set.seed(next_seed)

    for (disp in c(1, Inf)) {
      if (verbose) cat("ITER: ", iter, "DISPERSION: ", disp, "\n")

      sim_dat <- sim_data(ns, p, disp)
      dat <- sim_dat[[which_dat]]

      gb_res <- run_greedy_backfit(dat, Kmax = Kmax)
      all_res <- all_res |>
        bind_rows(next_tib(next_seed, disp, min(ns), "gb", Kmax, gb_res, sim_dat))
    }
  }

  return(all_res)
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


make_cosdist_plot <- function(spnmf_res, ebnmf_res, disp, method_name, xlabel, ylabel, title) {
  pens <- spnmf_res |> select(method) |> distinct()

  all_res <- spnmf_res |> mutate(source = method_name) |>
    bind_rows(ebnmf_res |> mutate(source = "EBNMF") |> select(-method) |> cross_join(pens)) |>
    group_by(method, Kmax, shape, metric_type, source) |>
    summarize(metric_na = anyNA(metric_val),
              metric_val = mean(1 - ifelse(metric_na, 0, metric_val)),
              .groups = "drop")

  plot_df <- all_res |>
    filter(str_detect(metric_type, "cosine")) |>
    mutate(LorF = factor(str_extract(metric_type, "L|F"), levels = c("L", "F"))) |>
    mutate(Component = paste("Component", str_extract(metric_type, "[0-9]+"))) |>
    mutate(source = factor(source, levels = c(method_name, "EBNMF")))

  if (disp == "Poisson") {
    plot_df <- plot_df |> filter(is.infinite(shape))
  } else { # "Overdispersed"
    plot_df <- plot_df |> filter(!is.infinite(shape))
  }

  p <- ggplot(plot_df, aes(x = method, y = metric_val, color = source, linetype = source)) +
    geom_line() +
    geom_point(data = plot_df |> filter(metric_na), aes(shape = metric_na), shape = 17) +
    scale_y_log10() +
    scale_color_manual(values = c("darkred", "dodgerblue")) +
    scale_linetype_manual(values = c("solid", "dashed")) +
    facet_grid(rows = vars(LorF), cols = vars(Component)) +
    labs(x = xlabel, y = ylabel, color = "Method") +
    # labs(caption = "Triangles indicate that components were removed in at least one trial.") +
    ggtitle(title) +
    guides(linetype = "none") +
    theme_bw() +
    theme(strip.text.y = element_text(angle = 0))

  return(p)
}

make_SNR_plot <- function(spnmf_res, ebnmf_res, disp, method_name, xlabel, ylabel) {
  pens <- spnmf_res |> select(method) |> distinct()

  all_res <- spnmf_res |> mutate(source = method_name) |>
    bind_rows(ebnmf_res |> mutate(source = "EBNMF") |> select(-method) |> cross_join(pens)) |>
    filter(str_detect(metric_type, "Scale")) |>
    group_by(method, Kmax, shape, metric_type, source, seed) |>
    mutate(metric_val = coalesce(metric_val, 0)) |>
    mutate(SNR = metric_val / (1 - sum(metric_val))) |>
    group_by(method, Kmax, shape, metric_type, source) |>
    summarize(SNR = mean(SNR), .groups = "drop") |>
    filter(SNR > 0)

  plot_df <- all_res |>
    mutate(Component = paste("Comp.", str_extract(metric_type, "[0-9]+"))) |>
    mutate(source = factor(source, levels = c(method_name, "EBNMF"))) |>
    filter(!str_detect(metric_type, "[7-9]"))

  if (disp == "Poisson") {
    plot_df <- plot_df |> filter(is.infinite(shape))
  } else { # "Overdispersed"
    plot_df <- plot_df |> filter(!is.infinite(shape))
  }

  p <- ggplot(plot_df, aes(x = method, y = SNR, color = source, linetype = source)) +
    geom_line() +
    scale_y_log10() +
    scale_color_manual(values = c("darkred", "dodgerblue")) +
    scale_linetype_manual(values = c("solid", "dashed")) +
    facet_grid(rows = vars(Component)) +
    labs(x = xlabel, y = ylabel, color = "Method") +
    guides(linetype = "none") +
    ggtitle("Noise-to-Signal") +
    theme_bw()

  return(p)
}

xlab <- "L1 penalty"
ylab <- "Cos. dist. from true L or F"
res1 <- ss1_rcppML_K4_res
res2 <- ss1_ebnmf_K4_res
p1 <- make_cosdist_plot(res1, res2, "Poisson", "RcppML", xlab, ylab,
                        "Poisson noise") +
  scale_x_sqrt()
p2 <- make_cosdist_plot(res1, res2, "Overdispersed", "RcppML", xlab, ylab,
                        "Overdispersed noise") +
  scale_x_sqrt()
plot_grid(p1, p2, nrow = 2)
ggsave("output/plots/ss1_rcppML_K4.pdf", width = 10, height = 6)

res1 <- ss1_nnlm_K4_res
p1 <- make_cosdist_plot(res1, res2, "Poisson", "NNLM", xlab, ylab,
                        "Poisson noise") +
  scale_x_sqrt()
p2 <- make_cosdist_plot(res1, res2, "Overdispersed", "NNLM", xlab, ylab,
                        "Overdispersed noise") +
  scale_x_sqrt()
plot_grid(p1, p2, nrow = 2)
ggsave("output/plots/ss1_nnlm_K4.pdf", width = 10, height = 6)

res1 <- ss1_hoyer_K4_res
xlab <- "Sparseness penalty"
p1 <- make_cosdist_plot(res1, res2, "Poisson", "Hoyer", xlab, ylab,
                        "Poisson noise")
p2 <- make_cosdist_plot(res1, res2, "Overdispersed", "Hoyer", xlab, ylab,
                        "Overdispersed noise")
plot_grid(p1, p2, nrow = 2)
ggsave("output/plots/ss1_hoyer_K4.pdf", width = 10, height = 6)

xlab <- "L1 penalty"
res1 <- ss1_rcppML_K6_res
res2 <- ss1_ebnmf_K6_res
p1a <- make_cosdist_plot(res1, res2, "Poisson", "RcppML", xlab, ylab,
                        "Poisson noise") +
  scale_x_sqrt() +
  guides(color = "none")
p1b <- make_SNR_plot(res1, res2, "Poisson", "RcppML", xlab, ylab)
p1 <- plot_grid(p1a, p1b, nrow = 1, rel_widths = c(0.65, 0.35))
p2a <- make_cosdist_plot(res1, res2, "Overdispersed", "RcppML", xlab, ylab,
                         "Overdispersed noise") +
  scale_x_sqrt() +
  guides(color = "none")
p2b <- make_SNR_plot(res1, res2, "Overdispersed", "RcppML", xlab, ylab)
p2 <- plot_grid(p2a, p2b, nrow = 1, rel_widths = c(0.65, 0.35))
plot_grid(p1, p2, nrow = 2)
ggsave("output/plots/ss1_rcppML_K6.pdf", width = 12, height = 8)

res1 <- ss1_nnlm_K6_res
p1a <- make_cosdist_plot(res1, res2, "Poisson", "NNLM", xlab, ylab,
                         "Poisson noise") +
  scale_x_sqrt() +
  guides(color = "none")
p1b <- make_SNR_plot(res1, res2, "Poisson", "NNLM", xlab, ylab)
p1 <- plot_grid(p1a, p1b, nrow = 1, rel_widths = c(0.65, 0.35))
p2a <- make_cosdist_plot(res1, res2, "Overdispersed", "NNLM", xlab, ylab,
                         "Overdispersed noise") +
  scale_x_sqrt() +
  guides(color = "none")
p2b <- make_SNR_plot(res1, res2, "Overdispersed", "NNLM", xlab, ylab)
p2 <- plot_grid(p2a, p2b, nrow = 1, rel_widths = c(0.65, 0.35))
plot_grid(p1, p2, nrow = 2)
ggsave("output/plots/ss1_nnlm_K6.pdf", width = 12, height = 8)

res1 <- ss1_hoyer_K6_res
xlab <- "Sparseness penalty"
p1a <- make_cosdist_plot(res1, res2, "Poisson", "Hoyer", xlab, ylab,
                         "Poisson noise") +
  guides(color = "none")
p1b <- make_SNR_plot(res1, res2, "Poisson", "Hoyer", xlab, ylab)
p1 <- plot_grid(p1a, p1b, nrow = 1, rel_widths = c(0.65, 0.35))
p2a <- make_cosdist_plot(res1, res2, "Overdispersed", "Hoyer", xlab, ylab,
                         "Overdispersed noise") +
  guides(color = "none")
p2b <- make_SNR_plot(res1, res2, "Overdispersed", "Hoyer", xlab, ylab)
p2 <- plot_grid(p2a, p2b, nrow = 1, rel_widths = c(0.65, 0.35))
plot_grid(p1, p2, nrow = 2)
ggsave("output/plots/ss1_Hoyer_K6.pdf", width = 12, height = 8)
