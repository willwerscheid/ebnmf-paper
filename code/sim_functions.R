## Fitting functions

run_nmf <- function(Y, k, ntrials = 10) {
  cat("Running NMF")
  t <- system.time({
    best_mse <- Inf
    for (i in 1:ntrials) {
      cat(".")
      set.seed(i)
      next_res <- NNLM::nnmf(Y, k = k, verbose = 0)
      if (min(next_res$mse) < best_mse) {
        best_mse <- min(next_res$mse)
        best_res <- next_res
      }
    }
  })
  cat("\n")

  return(list(t = t, fit = best_res))
}

run_ebnmf_from_nmf <- function(Y, nmf_res, var_type, maxiter = 2000) {
  cat("Running EBNMF from NMF:\n")
  t <- system.time({
    fl <- flash_init(Y, var_type = var_type) |>
      flash_factors_init(list(nmf_res$W, t(nmf_res$H)), ebnm_fn = ebnm_point_exponential) |>
      flash_backfit(maxiter = maxiter, verbose = 0) |>
      flash_nullcheck(verbose = 0)
  })

  return(list(t = t, fit = fl))
}

run_greedy_backfit <- function(Y, Kmax, var_type = 2) {
  cat("Running greedy-backfit EBNMF")
  t <- system.time({
    fl <- flash(Y, var_type = var_type, greedy_Kmax = Kmax, ebnm_fn = ebnm_point_exponential) |>
      flash_backfit() |>
      flash_nullcheck()
    keep_going <- fl$n_factors < Kmax
    while (keep_going) {
      cat(".")
      current_n <- fl$n_factors
      fl <- fl |>
        flash_greedy(Kmax = Kmax - current_n, ebnm_fn = ebnm_point_exponential) |>
        flash_backfit() |>
        flash_nullcheck()
      if (fl$n_factors == current_n | fl$n_factors == Kmax) {
        keep_going <- FALSE
      }
    }
    cat("\n")
  })

  return(list(t = t, fit = fl))
}

run_alternating <- function(Y, Kmax, var_type = 2) {
  cat("Running alternating EBNMF")
  t <- system.time({
    fl <- flash_init(Y, var_type = var_type) |>
      flash_set_verbose(0)
    keep_going <- TRUE
    while(keep_going) {
      cat(".")
      current_n <- fl$n_factors
      fl <- fl |>
        flash_greedy(ebnm_fn = ebnm_point_exponential) |>
        flash_backfit(maxiter = 10)
      if (fl$n_factors == current_n | fl$n_factors == Kmax) {
        keep_going <- FALSE
      }
    }
    cat("\n")
    fl <- fl |>
      flash_backfit(maxiter = 2000, verbose = 0) |>
      flash_nullcheck(verbose = 0)
  })

  return(list(t = t, fit = fl))
}

## Evaluation metrics

calc_metrics <- function(res, sim_dat) {
  fit <- res$fit
  if (inherits(fit, "flash")) {
    LDF <- ldf(fit, type = "f")
    LL <- t(t(LDF$L) * LDF$D)
    FF <- t(t(LDF$F) * LDF$D)
  } else {
    Wscale <- sqrt(apply(fit$W, 2, function(x) sum(x^2)))
    Hscale <- sqrt(apply(fit$H, 1, function(x) sum(x^2)))
    D <- Wscale * Hscale
    LL <- t(t(fit$W) / sqrt(Wscale) * sqrt(Hscale))
    FF <- t(fit$H / sqrt(Hscale) * sqrt(Wscale))
  }

  LL_cors <- FF_cors <- rep(NA, trueK)
  used_cols <- numeric(0)
  LL_cormat <- cor(sim_dat$L, LL)
  FF_cormat <- cor(sim_dat$F, FF)
  for (i in 1:min(trueK, ncol(LL))) {
    rowmax <- which.max(apply(abs(LL_cormat), 1, max))
    colmax <- which.max(apply(abs(LL_cormat), 2, max))
    LL_cors[rowmax] <- LL_cormat[rowmax, colmax]
    FF_cors[rowmax] <- FF_cormat[rowmax, colmax]
    LL_cormat[, colmax] <- 0
    LL_cormat[rowmax, ] <- 0
    used_cols <- c(used_cols, colmax)
  }

  L_scales <- apply(LL, 2, function(x) sum(x^2))
  L_scales <- L_scales / sum(L_scales)
  unmatched_scales <- L_scales[-used_cols]
  unmatched_scales <- sort(unmatched_scales, decreasing = TRUE)
  unmatched_scales <- c(unmatched_scales, rep(NA, highK - trueK - length(unmatched_scales)))

  all_metrics <- c(LL_cors, FF_cors, unmatched_scales)
  names(all_metrics) <- c(
    paste0("LLcor", 1:trueK),
    paste0("FFcor", 1:trueK),
    paste0("Scale", (trueK + 1):highK)
  )

  return(all_metrics)
}

next_tib <- function(shape, ns, method, Kmax, res, sim_dat) {
  metrics <- calc_metrics(res, sim_dat)
  return(tibble(
    method = method,
    Kmax = Kmax,
    varied_n = ns[4],
    shape = shape,
    metric_type = c("t_elapsed", names(metrics)),
    metric_val = c(res$t[3], metrics)
  ))
}

## Plotting functions

make_corplot <- function(res, fill_label) {
  plot_df <- res |>
    mutate(varied_n = factor(varied_n),
           shape = factor(round(shape, 2)),
           method = factor(method, levels = lvls)) |>
    mutate(metric_val = pmax(pmin(metric_val, 0.99999), 0.1))
  p <- ggplot(plot_df, aes(x = varied_n, y = shape, fill = metric_val)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "blue", na.value = "black", transform = "logit",
                        breaks = c(0.5, 0.9, 0.99, 0.999, 0.9999)) +
    facet_grid(rows = vars(metric_type), cols = vars(method), scales = "free_x") +
    labs(x = xlabel, y = "Shape of gamma prior on factors", fill = fill_label) +
    theme(axis.text.x = element_text(angle = 45))
  return(p)
}

make_scaleplot <- function(res, fill_label) {
  plot_df <- res |>
    mutate(varied_n = factor(varied_n),
           shape = factor(round(shape, 2)),
           method = factor(method, levels = lvls)) |>
    mutate(metric_val = pmax(metric_val, 0.00001))
  p <- ggplot(plot_df, aes(x = varied_n, y = shape, fill = metric_val)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "red", na.value = "white", transform = "logit",
                        breaks = c(.0001, .001, .01, .1)) +
    facet_grid(rows = vars(metric_type), cols = vars(method), scales = "free_x") +
    labs(x = xlabel, y = "Shape of gamma prior on factors", fill = "Scale") +
    theme(axis.text.x = element_text(angle = 45))
  return(p)
}

align_cols <- function(plotmat, refmat) {
  if (ncol(refmat) > ncol(plotmat)) {
    refmat <- refmat[, 1:ncol(plotmat)]
  }
  idx <- rep(NA, ncol(plotmat))
  cor_mat <- cor(plotmat, refmat)
  cor_which <- apply(cor_mat, 1, which.max)
  all_k <- rev(order(apply(cor_mat, 1, max)))
  for (k in all_k) {
    next_idx <- cor_which[k]
    if (is.na(idx[next_idx])) {
      idx[next_idx] <- k
    }
  }
  idx[which(is.na(idx))] <- setdiff(1:ncol(plotmat), idx)
  return(idx)
}

make_structure_plot <- function(LL, kset, title) {
  structure_plot(LL[, kset], grouping = rep(LETTERS[1:length(ns)], times = ns), gap = 20,
                 topics = rev(1:ncol(LL)), loadings_order = 1:sum(ns),
                 colors = RColorBrewer::brewer.pal(12, "Set3")) +
    labs(y = "") +
    ggtitle(title) +
    guides(fill = "none", color = "none") +
    theme(axis.text.x = element_blank()) +
    theme(axis.text.y = element_blank()) +
    theme(plot.title.position = "plot")
}

plot_fl <- function(fl, refmat, title) {
  LDF <- ldf(fl, type = "f")
  LL <- t(t(LDF$L) * LDF$D)
  kset <- align_cols(LL, refmat)
  make_structure_plot(LL, kset, title)
}

plot_nmf <- function(res, refmat, title) {
  Wscale <- sqrt(apply(res$W, 2, function(x) sum(x^2)))
  Hscale <- sqrt(apply(res$H, 1, function(x) sum(x^2)))
  D <- Wscale * Hscale
  LL <- t(t(res$W) / sqrt(Wscale) * sqrt(Hscale))
  kset <- align_cols(LL, refmat)
  make_structure_plot(LL, kset, title)
}
