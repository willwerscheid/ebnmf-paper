run_RcppML_sims <- function(sim_data, which_dat, Kmax, niter = 5, verbose = FALSE) {
  all_res <- tibble()
  next_seed <- 0
  for (iter in 1:niter) {
    if (verbose) cat("ITER: ", iter, "\n")

    next_seed <- next_seed + 1
    set.seed(next_seed)

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

  return(all_res)
}

run_KimPark_sims <- function(matlab, which_dat, Kmax, L1pens, L2pen = 1, niter = 5, verbose = FALSE) {
  all_res <- tibble()
  next_seed <- 0
  for (iter in 1:niter) {
    next_seed <- next_seed + 1
    set.seed(next_seed)

    for (disp in c(1, Inf)) {
      if (verbose) cat("RARE N: ", 10, "DISPERSION: ", disp, "\n")

      sim_dat <- sim_data(ns, p, disp)
      dat <- sim_dat[[which_dat]]

      for (pen in l1pens) {
        if (verbose) cat(" penalty = ", pen, "...\n")

        spnmf_res <- run_KimPark_sparse_nmf(matlab, dat, Kmax, pen, L2pen, 1:niter, verbose)
        all_res <- all_res |>
          bind_rows(next_tib(next_seed, disp, min(ns), pen, Kmax, spnmf_res, sim_dat))
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
          all_t <- all_t + max(hoyer_res$t, na.rm = TRUE)
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

make_SNR_plot <- function(spnmf_res, ebnmf_res, disp, method_name, xlabel, ylabel, trunc_at) {
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

  if (trunc_at < 10) {
    complab <- "Component"
  } else {
    complab <- "Comp."
  }
  plot_df <- all_res |>
    mutate(Component = paste(complab, str_extract(metric_type, "[0-9]+"))) |>
    mutate(source = factor(source, levels = c(method_name, "EBNMF"))) |>
    mutate(Component = factor(Component, levels = paste(complab, 1:trunc_at))) |>
    filter(!is.na(Component))

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
    facet_grid(cols = vars(Component)) +
    labs(x = xlabel, y = ylabel, color = "Method") +
    guides(linetype = "none") +
    theme_bw()

  return(p)
}

make_trueK_plot <- function(res1, res2, method_name, xlab, ylab,
                            scale_sqrt = FALSE, brks = NULL) {
  p1 <- make_cosdist_plot(
    res1, res2, "Poisson", method_name, xlab, "Cos. dist. from true L or F", "Poisson noise"
  )
  p2 <- make_cosdist_plot(
    res1, res2, "Overdispersed", method_name, xlab, "Cos. dist. from true L or F", "Overdispersed noise"
  )
  if (scale_sqrt) {
    p1 <- p1 + scale_x_sqrt(breaks = brks)
    p2 <- p2 + scale_x_sqrt(breaks = brks)
  } else if (!is.null(brks)) {
    p1 <- p1 + scale_x_continuous(breaks = brks)
    p2 <- p2 + scale_x_continuous(breaks = brks)
  }
  plot_grid(p1, p2, nrow = 2)
}

make_overK_plot <- function(res1, res2, method_name, xlab, ylab, rel_widths, trunc_at,
                            scale_sqrt = FALSE, brks = NULL) {
  p1a <- make_cosdist_plot(
    res1, res2, "Poisson", method_name, xlab, "Cos. dist. from true L or F", "Poisson noise"
  ) +
    guides(color = "none")
  p1b <- make_SNR_plot(res1, res2, "Poisson", method_name, xlab, "Noise-to-signal ratio", trunc_at)
  p2a <- make_cosdist_plot(
    res1, res2, "Overdispersed", method_name, xlab, "Cos. dist. from true L or F", "Overdispersed noise") +
    guides(color = "none")
  p2b <- make_SNR_plot(res1, res2, "Overdispersed", method_name, xlab, "Noise-to-signal ratio", trunc_at)
  if (scale_sqrt) {
    p1a <- p1a + scale_x_sqrt(breaks = brks)
    p1b <- p1b + scale_x_sqrt(breaks = brks)
    p2a <- p2a + scale_x_sqrt(breaks = brks)
    p2b <- p2b + scale_x_sqrt(breaks = brks)
  } else if (!is.null(brks)) {
    p1a <- p1a + scale_x_continuous(breaks = brks)
    p1b <- p1b + scale_x_continuous(breaks = brks)
    p2a <- p2a + scale_x_continuous(breaks = brks)
    p2b <- p2b + scale_x_continuous(breaks = brks)
  }
  p1b_wmargin <- suppressWarnings({
    plot_grid("", p1b, "", ncol = 1, rel_heights = c(0.2, 0.6, 0.2))
  })
  p2b_wmargin <- suppressWarnings({
    plot_grid("", p2b, "", ncol = 1, rel_heights = c(0.2, 0.6, 0.2))
  })
  p1 <- plot_grid(p1a, p1b_wmargin, rel_widths = rel_widths)
  p2 <- plot_grid(p2a, p2b_wmargin, rel_widths = rel_widths)
  plot_grid(p1, p2, nrow = 2)
}
