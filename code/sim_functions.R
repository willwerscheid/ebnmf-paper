## Simulation functions

sim_F <- function(p, k, gamma_shape, gamma_scale, n_anchor_words) {
  F <- matrix(rgamma(p * k, shape = gamma_shape, scale = gamma_scale), nrow = p, ncol = k)

  # Anchor words
  for (i in 1:k) {
    F[((i - 1) * n_anchor_words + 1):(i * n_anchor_words), i] <- 1
    F[((i - 1) * n_anchor_words + 1):(i * n_anchor_words), setdiff(1:k, i)] <- 0
  }

  return(F)
}

sim_X <- function(L, F, dispersion = Inf) {
  mu <- L %*% t(F)
  if (is.infinite(dispersion)) {
    X <- matrix(rpois(length(mu), expm1(mu)), nrow = nrow(mu), ncol = ncol(mu))
  } else {
    X <- matrix(rnbinom(length(mu), mu = expm1(mu), size = dispersion), nrow = nrow(mu), ncol = ncol(mu))
  }
  return(X)
}

return_sim_data <- function(X, L, F, pops) {
  # Make sure there aren't any all-zero columns:
  F <- F[apply(X, 2, sum) > 0, ]
  X <- X[, apply(X, 2, sum) > 0]

  Xnorm <- t(t(X) / colSums(X)) * median(colSums(X))
  Xlibnorm <- X / rowSums(X) * median(rowSums(X))

  Y <- log1p(X)
  Ynorm <- log1p(Xnorm)
  Ylibnorm <- log1p(Xlibnorm)

  rownames(Y) <- paste0("sample", 1:nrow(Y))
  colnames(Y) <- paste0("feature", 1:ncol(Y))

  return(list(Y = Y, Ynorm = Ynorm, Ylibnorm = Ylibnorm, L = L, F = F, pops = pops))
}

## Fitting functions

run_nmf <- function(Y, k, seeds = 1:10, method = "scd", verbose = FALSE) {
  if (verbose) cat("Running NMF")
  t <- system.time({
    all_mse <- numeric(length(seeds))
    best_mse <- Inf
    # The seeds parameter controls the number of trials:
    for (i in 1:length(seeds)) {
      if (verbose) cat(".")
      set.seed(seeds[i])
      next_res <- NNLM::nnmf(Y, k = k, verbose = 0, method = method)
      all_mse[i] <- min(next_res$mse)
      if (min(next_res$mse) < best_mse) {
        best_mse <- min(next_res$mse)
        best_res <- next_res
      }
    }
  })
  if (verbose) cat("\n")

  return(list(t = t, fit = best_res, all_mse = all_mse))
}

run_sparse_nmf <- function(Y, k, L1pen, seeds = 1, verbose = FALSE) {
  if (verbose) cat("Running Sparse NMF, L1 penalty =", L1pen)
  t <- system.time({
    all_mse <- numeric(length(seeds))
    best_mse <- Inf
    # The seeds parameter controls the number of trials:
    for (i in 1:length(seeds)) {
      if (verbose) cat(".")
      set.seed(seeds[i])
      next_res <- NNLM::nnmf(Y, k = k, alpha = c(0, 0, L1pen), beta = c(0, 0, L1pen), verbose = 0)
      all_mse[i] <- min(next_res$mse)
      if (min(next_res$mse) < best_mse) {
        best_mse <- min(next_res$mse)
        best_res <- next_res
      }
    }
  })
  if (verbose) cat("\n")

  return(list(t = t, fit = best_res, all_mse = all_mse))
}

run_RcppML_sparse_nmf <- function(Y, k, L1pen, seeds = 1, verbose = FALSE) {
  if (verbose) cat("Running Sparse NMF, L1 penalty =", L1pen)
  t <- system.time({
    all_mse <- numeric(length(seeds))
    best_mse <- Inf
    # The seeds parameter controls the number of trials:
    for (i in 1:length(seeds)) {
      if (verbose) cat(".")
      next_res <- RcppML::nmf(Y, k = k, L1 = c(L1pen, 0), verbose = FALSE, seed = seeds[i])
      all_mse[i] <- mean((Y - next_res$w %*% (next_res$h * next_res$d))^2)
      if (all_mse[i] < best_mse) {
        best_mse <- all_mse[i]
        best_res <- next_res
      }
    }
  })

  out_res <- list()
  out_res$W <- t(t(best_res$w) * sqrt(best_res$d))
  out_res$H <- best_res$h * sqrt(best_res$d)

  out_res$W <- out_res$W[, colSums(best_res$w) > 0]
  out_res$H <- out_res$H[colSums(best_res$w) > 0, ]

  return(list(t = t, fit = out_res, all_mse = all_mse))
}

run_ebnmf_from_nmf <- function(Y, nmf_res, var_type = 2, ebnm_fn = ebnm_point_exponential,
                               maxiter = 2000, verbose = FALSE) {
  if (verbose) cat("Running EBNMF from NMF.\n")
  t <- system.time({
    fl <- flash_init(Y, var_type = var_type) |>
      flash_factors_init(list(nmf_res$W, t(nmf_res$H)), ebnm_fn = ebnm_fn) |>
      flash_backfit(maxiter = maxiter, verbose = 0) |>
      flash_nullcheck(verbose = 0)
  })

  return(list(t = t, fit = fl))
}

run_greedy_backfit <- function(Y, Kmax, var_type = 2, ebnm_fn = ebnm_point_exponential,
                               alternating = TRUE, verbose = FALSE) {
  if (verbose) cat("Running greedy-backfit EBNMF")
  t <- system.time({
    fl <- flash(Y, var_type = var_type, greedy_Kmax = Kmax, ebnm_fn = ebnm_fn, verbose = 0) |>
      flash_backfit(verbose = 0) |>
      flash_nullcheck(verbose = 0)
    keep_going <- fl$n_factors < Kmax
    while (keep_going) {
      if (verbose) cat(".")
      current_n <- fl$n_factors
      fl <- fl |>
        flash_greedy(Kmax = Kmax - current_n, ebnm_fn = ebnm_fn, verbose = 0) |>
        flash_backfit(verbose = 0) |>
        flash_nullcheck(verbose = 0)
      if (!alternating | fl$n_factors == current_n | fl$n_factors == Kmax) {
        keep_going <- FALSE
      }
    }
    if (verbose) cat("\n")
  })

  return(list(t = t, fit = fl))
}

run_alternating <- function(Y, Kmax, var_type = 2, ebnm_fn = ebnm_point_exponential, verbose = FALSE) {
  if (verbose) cat("Running alternating EBNMF")
  t <- system.time({
    fl <- flash_init(Y, var_type = var_type) |>
      flash_set_verbose()
    keep_going <- TRUE
    while(keep_going) {
      if (verbose) cat(".")
      current_n <- fl$n_factors
      fl <- fl |>
        flash_greedy(ebnm_fn = ebnm_fn, verbose = 0) |>
        flash_backfit(maxiter = 10, verbose = 0)
      if (fl$n_factors == current_n | fl$n_factors == Kmax) {
        keep_going <- FALSE
      }
    }
    if (verbose) cat("\n")
    fl <- fl |>
      flash_backfit(maxiter = 2000, verbose = 0) |>
      flash_nullcheck(verbose = 0)
  })

  return(list(t = t, fit = fl))
}

## CV style approaches

run_RcppML_cv <- function(Y, k, L1pens, nfolds = 100, ntrials = 10, verbose = FALSE) {
  folds <- sample(1:nfolds, length(Y), replace = TRUE)
  L1pen_mse <- numeric(length(L1pens))
  for (i in 1:length(L1pens)) {
    if (verbose) cat("L1 penalty =", L1pens[i], "\n")
    all_mse <- numeric(ntrials)
    for (fold in 1:ntrials) {
      if (verbose) cat(". Fold =", fold, "\n")
      dat <- Y + 1e-4
      dat[folds == fold] <- 0
      dat <- Matrix(dat, sparse = TRUE)
      nmf_res <- RcppML::nmf(dat, k = k, L1 = c(L1pens[i], 0), mask_zeros = TRUE, verbose = FALSE)
      imp_dat <- nmf_res$w %*% diag(nmf_res$d) %*% nmf_res$h
      mse <- mean((imp_dat[folds == fold] - Y[folds == fold])^2)
      all_mse[fold] <- mse
    }
    if (verbose) cat("MSE =", mean(all_mse), "\n")
    L1pen_mse[i] <- mean(all_mse)
  }
  return(tibble(L1pen = L1pens, RMSE = sqrt(L1pen_mse)))
}

run_ebnmf_cv <- function(Y, k, nfolds = 100, ntrials = 10, verbose = FALSE) {
  folds <- sample(1:nfolds, length(Y), replace = TRUE)
  all_mse <- numeric(ntrials)
  for (fold in 1:ntrials) {
    if (verbose) cat("FOLD =", fold, "\n")
    dat <- Y
    dat[folds == fold] <- NA
    ebnmf_res <- run_alternating(dat, Kmax = k, var_type = 2, verbose = verbose)
    imp_dat <- fitted(ebnmf_res$fit)
    mse <- mean((imp_dat[folds == fold] - Y[folds == fold])^2)
    all_mse[fold] <- mse
  }
  if (verbose) cat("MSE =", mean(all_mse), "\n")
  return(sqrt(mean(all_mse)))
}

## Evaluation metrics

calc_metrics <- function(res, sim_dat) {
  Lscale <- sqrt(apply(sim_dat$L, 2, function(x) sum(x^2)))
  Fscale <- sqrt(apply(sim_dat$F, 2, function(x) sum(x^2)))
  true_L <- t(t(sim_dat$L) / Lscale)
  true_F <- t(t(sim_dat$F) / Fscale)

  LL_cosine <- FF_cosine <- rep(NA, ncol(sim_dat$L))

  fit <- res$fit
  if (inherits(fit, "flash")) {
    LDF <- ldf(fit, type = "f")
    LL <- LDF$L
    FF <- LDF$F
    D <- LDF$D
  } else if (is.matrix(fit$W)) {
    Wscale <- sqrt(apply(fit$W, 2, function(x) sum(x^2)))
    Hscale <- sqrt(apply(fit$H, 1, function(x) sum(x^2)))
    LL <- t(t(fit$W) / Wscale)
    FF <- t(fit$H / Hscale)
    D <- Wscale * Hscale
  } else {
    return(c(LL_cosine, FF_cosine))
  }

  if (ncol(LL) == 0) {
    return(c(LL_cosine, FF_cosine))
  }

  used_cols <- numeric(0)
  LL_cosmat <- crossprod(true_L, LL)
  FF_cosmat <- crossprod(true_F, FF)
  for (i in 1:min(ncol(true_L), ncol(LL))) {
    rowmax <- which.max(apply(abs(LL_cosmat), 1, max))
    colmax <- which.max(apply(abs(LL_cosmat), 2, max))
    LL_cosine[rowmax] <- LL_cosmat[rowmax, colmax]
    FF_cosine[rowmax] <- FF_cosmat[rowmax, colmax]
    LL_cosmat[, colmax] <- 0
    LL_cosmat[rowmax, ] <- 0
    used_cols <- c(used_cols, colmax)
  }

  unmatched_scales <- D[-used_cols] / sum(D)
  unmatched_scales <- sort(unmatched_scales, decreasing = TRUE)

  names(LL_cosine) <- paste0("LLcosine", 1:length(LL_cosine))
  names(FF_cosine) <- paste0("FFcosine", 1:length(FF_cosine))
  if (length(unmatched_scales) > 0) {
    names(unmatched_scales) <- paste0("Scale", (ncol(sim_dat$L) + 1):(ncol(sim_dat$L) + length(unmatched_scales)))
  }

  all_metrics <- c(LL_cosine, FF_cosine, unmatched_scales)
  return(all_metrics)
}

next_tib <- function(seed, shape, varied_n, method, Kmax, res, sim_dat) {
  metrics <- calc_metrics(res, sim_dat)
  return(tibble(
    seed = seed,
    method = method,
    Kmax = Kmax,
    varied_n = varied_n,
    shape = shape,
    metric_type = c("t_elapsed", names(metrics)),
    metric_val = c(res$t[3], metrics)
  ))
}

## Plotting functions

make_cosplot <- function(res, metric, xlabel, ylabel, fill_label, cutoff = 0.99999) {
  res <- res |>
    filter(str_starts(metric_type, metric))
  plot_df <- res |>
    mutate(varied_n = factor(varied_n),
           shape = factor(round(shape, 2))) |>
    mutate(metric_val = pmax(pmin(metric_val, cutoff), 0.1))
  p <- ggplot(plot_df, aes(x = varied_n, y = shape, fill = metric_val)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "blue", na.value = "gray", transform = "logit",
                        breaks = c(0.5, 1 - 10^seq(-1, ceiling(log10(1 - cutoff)), by = -1)),
                        limits = c(0.1, cutoff)) +
    facet_grid(rows = vars(metric_type), cols = vars(method), scales = "free_x") +
    labs(x = xlabel, y = ylabel, fill = fill_label) +
    theme(axis.text.x = element_text(angle = 45, size = 6)) +
    theme(strip.text.x = element_text(size = 4), strip.text.y = element_text(size = 6))
  return(p)
}

make_scaleplot <- function(res, xlabel, ylabel, fill_label, cutoff = 1000) {
  res <- res |>
    filter(str_starts(metric_type, "Scale")) |>
    mutate(metric_type = factor(metric_type, levels = paste0("Scale", 1:12))) |>
    group_by(varied_n, shape, method) |>
    mutate(SNR = (1 - sum(metric_val)) / metric_val)
  plot_df <- res |>
    mutate(varied_n = factor(varied_n),
           shape = factor(round(shape, 2))) |>
    mutate(SNR = pmin(SNR, cutoff))
  p <- ggplot(plot_df, aes(x = varied_n, y = shape, fill = SNR)) +
    geom_tile() +
    scale_fill_gradient(low = "red", high = "white", na.value = "white", transform = "log",
                        breaks = 10^(seq(0, floor(log10(cutoff)), by = 1)),
                        limits = c(1, cutoff)) +
    facet_grid(rows = vars(metric_type), cols = vars(method), scales = "free_x") +
    labs(x = xlabel, y = ylabel, fill = fill_label) +
    theme(axis.text.x = element_text(angle = 45, size = 6)) +
    theme(strip.text.x = element_text(size = 4), strip.text.y = element_text(size = 6))
  return(p)
}

make_timingplot <- function(res, xlabel, ylabel) {
  plot_df <- res |> filter(str_starts(metric_type, "t_elapsed")) |>
    mutate(varied_n = factor(varied_n),
           shape = factor(round(shape, 2)),
           Kmax = factor(paste0("Kmax = ", Kmax))) |>
    group_by(varied_n, shape, Kmax, method) |>
    summarize(metric_val = mean(metric_val, na.rm = TRUE))
  ggplot(plot_df, aes(x = varied_n, y = shape, fill = metric_val)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "darkgreen", na.value = "black", transform = "log10") +
    facet_grid(rows = vars(Kmax), cols = vars(method), scales = "free_x") +
    labs(x = xlabel, y = ylabel, fill = "Time elapsed (s)") +
    theme(axis.text.x = element_text(angle = 45, size = 6)) +
    theme(strip.text.x = element_text(size = 4), strip.text.y = element_text(size = 6))
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

make_structure_plot <- function(LL, pops, kset, title) {
  structure_plot(LL[, kset], grouping = pops, gap = 20,
                 topics = rev(1:ncol(LL)), loadings_order = 1:nrow(LL),
                 colors = RColorBrewer::brewer.pal(12, "Set3")) +
    labs(y = "") +
    ggtitle(title) +
    guides(fill = "none", color = "none") +
    theme(axis.text.x = element_blank()) +
    theme(axis.text.y = element_blank()) +
    theme(plot.title.position = "plot")
}

plot_fl <- function(fl, refmat, pops, title) {
  LDF <- ldf(fl, type = "f")
  LL <- t(t(LDF$L) * LDF$D)
  kset <- align_cols(LL, refmat)
  make_structure_plot(LL, pops, kset, title)
}

plot_nmf <- function(fit, refmat, pops, title) {
  Wscale <- sqrt(apply(fit$W, 2, function(x) sum(x^2)))
  Hscale <- sqrt(apply(fit$H, 1, function(x) sum(x^2)))
  D <- Wscale * Hscale
  LL <- t(t(fit$W) / sqrt(Wscale) * sqrt(Hscale))
  kset <- align_cols(LL, refmat)
  make_structure_plot(LL, pops, kset, title)
}
