source("./code/sim_functions.R")
options(matlab = "/Applications/MATLAB_R2025a.app/bin/matlab")


sim_data <- function(ns, p, colwise_noise, component_scales = c(4, 2, 4, 1, 4)) {
  pops <- rep(LETTERS[1:length(ns)], times = ns)

  L <- matrix(0, nrow = sum(ns), ncol = 5)
  L[, 1] <- 1 # "background" factor
  L[pops == "A", 2] <- 1
  L[pops == "B", 3] <- 1
  L[pops == "C", 4] <- 1
  L[pops == "D", 5] <- 1

  L <- vary_library_sizes(L)
  F <- sim_F(p, ncol(L))
  F <- t(t(F) * component_scales)
  X <- sim_X(L, F, colwise_noise = colwise_noise)
  return_sim_data(X, L, F, pops)
}

ns <- c(1250, 125, 125, 20)
p <- 2000
set.seed(1)
example_sim <- sim_data(ns, p, Inf)

Y <- example_sim$Y
L <- example_sim$L
pops <- example_sim$pops


nmf_k4_res <- run_RcppML_sparse_nmf(Y, k = 4, L1pen = 0, seeds = 1:10)
nmf_k4_worst_res <- run_RcppML_sparse_nmf(Y, k = 4, L1pen = 0, seeds = which.max(nmf_k4_res$all_mse))
nmf_k4_p1 <- plot_nmf(
  nmf_k4_res$fit, L, pops,
  paste0("K = 4, best fit (MSE = ", round(min(nmf_k4_res$all_mse), 4), ")"),
  gap = 50
)
nmf_k4_p2 <- plot_nmf(
  nmf_k4_worst_res$fit, L, pops,
  paste0("K = 4, worst fit (MSE = ", round(nmf_k4_worst_res$all_mse, 4), ")"),
  gap = 50
)

nmf_k5_res <- run_RcppML_sparse_nmf(Y, k = 5, L1pen = 0, seeds = 1:10)
nmf_k5_worst_res <- run_RcppML_sparse_nmf(Y, k = 5, L1pen = 0, seeds = which.max(nmf_k5_res$all_mse))
nmf_k5_p1 <- plot_nmf(
  nmf_k5_res$fit, L, pops,
  paste0("K = 5, best fit (MSE = ", round(min(nmf_k5_res$all_mse), 4), ")"),
  gap = 50
)
nmf_k5_p2 <- plot_nmf(
  nmf_k5_worst_res$fit, L, pops,
  paste0("K = 5, worst fit (MSE = ", round(nmf_k5_worst_res$all_mse, 4), ")"),
  gap = 50
)

nmf_k8_res <- run_RcppML_sparse_nmf(Y, k = 8, L1pen = 0, seeds = 1:10)
nmf_k8_worst_res <- run_RcppML_sparse_nmf(Y, k = 8, L1pen = 0, seeds = which.max(nmf_k8_res$all_mse))
nmf_k8_p1 <- plot_nmf(
  nmf_k8_res$fit, L, pops,
  paste0("K = 8, best fit (MSE = ", round(min(nmf_k8_res$all_mse), 4), ")"),
  gap = 50
)
nmf_k8_p2 <- plot_nmf(
  nmf_k8_worst_res$fit, L, pops,
  paste0("K = 8, worst fit (MSE = ", round(nmf_k8_worst_res$all_mse, 4), ")"),
  gap = 50
)

plot_grid(nmf_k4_p1, nmf_k4_p2, nmf_k5_p1, nmf_k5_p2, nmf_k8_p1, nmf_k8_p2,
          nrow = 2, byrow = FALSE)
ggsave("output/plots/ex_nmf.pdf", width = 7, height = 5)


setwd("./matlab/")
Matlab$startServer(matlab = "/Applications/MATLAB_R2025a.app/bin/matlab")
matlab <- Matlab()
open(matlab)
setVerbose(matlab, threshold = 100000)

KimPark_k4_plots <- list()
KimPark_k4_t <- numeric()
for (L1pen in c(1, 10, 100, 1000)) {
  next_res <- run_KimPark_sparse_nmf(matlab, Y, k = 4, L1pen = L1pen, L2pen = 1, seeds = 1:10)
  if (!is.null(next_res$fit$W)) {
    next_p <- plot_nmf(next_res$fit, L, pops, paste0("K = 4, L1 penalty = ", L1pen), gap = 50)
    KimPark_k4_plots <- c(KimPark_k4_plots, list(next_p))
    KimPark_k4_t <- c(KimPark_k4_t, next_res$t[3])
  }
}

KimPark_k5_plots <- list()
KimPark_k5_t <- numeric()
for (L1pen in c(1, 10, 100, 1000)) {
  next_res <- run_KimPark_sparse_nmf(matlab, Y, k = 5, L1pen = L1pen, L2pen = 1, seeds = 1:10)
  if (!is.null(next_res$fit$W)) {
    next_p <- plot_nmf(next_res$fit, L, pops, paste0("K = 5, L1 penalty =  ", L1pen), gap = 50)
    KimPark_k5_plots <- c(KimPark_k5_plots, list(next_p))
    KimPark_k5_t <- c(KimPark_k5_t, next_res$t[3])
  }
}

KimPark_k8_plots <- list()
KimPark_k8_t <- numeric()
for (L1pen in c(1, 10, 100, 1000)) {
  next_res <- run_KimPark_sparse_nmf(matlab, Y, k = 8, L1pen = L1pen, L2pen = 1, seeds = 1:10)
  if (!is.null(next_res$fit$W)) {
    next_p <- plot_nmf(next_res$fit, L, pops, paste0("K = 8, L1 penalty =  ", L1pen), gap = 50)
    KimPark_k8_plots <- c(KimPark_k8_plots, list(next_p))
    KimPark_k8_t <- c(KimPark_k8_t, next_res$t[3])
  }
}

Hoyer_k4_plots <- list()
Hoyer_k4_t <- numeric()
for (L1pen in seq(0.1, 0.7, by = 0.2)) {
  next_res <- run_Hoyer_sparse_nmf(matlab, Y, k = 4, pen = L1pen, seeds = 1:10)
  if (!is.null(next_res$fit$W)) {
    next_p <- plot_nmf(next_res$fit, L, pops, paste0("K = 4, Sparseness = ", L1pen), gap = 50)
    Hoyer_k4_plots <- c(Hoyer_k4_plots, list(next_p))
    Hoyer_k4_t <- c(Hoyer_k4_t, next_res$t[3])
  }
}

Hoyer_k5_plots <- list()
Hoyer_k5_t <- numeric()
for (L1pen in seq(0.1, 0.7, by = 0.2)) {
  next_res <- run_Hoyer_sparse_nmf(matlab, Y, k = 5, pen = L1pen, seeds = 1:10)
  if (!is.null(next_res$fit$W)) {
    next_p <- plot_nmf(next_res$fit, L, pops, paste0("K = 5, Sparseness = ", L1pen), gap = 50)
    Hoyer_k5_plots <- c(Hoyer_k5_plots, list(next_p))
    Hoyer_k5_t <- c(Hoyer_k5_t, next_res$t[3])
  }
}

Hoyer_k8_plots <- list()
Hoyer_k8_t <- numeric()
for (L1pen in seq(0.1, 0.7, by = 0.2)) {
  next_res <- run_Hoyer_sparse_nmf(matlab, Y, k = 8, pen = L1pen, seeds = 1:10)
  if (!is.null(next_res$fit$W)) {
    next_p <- plot_nmf(next_res$fit, L, pops, paste0("K = 8, Sparseness = ", L1pen), gap = 50)
    Hoyer_k8_plots <- c(Hoyer_k8_plots, list(next_p))
    Hoyer_k8_t <- c(Hoyer_k8_t, next_res$t[3])
  }
}

close(matlab)
setwd("../")


nmf_k4 <- nmf_k4_p1 + ggtitle("K = 4, Vanilla NMF")
nmf_k5 <- nmf_k5_p1 + ggtitle("K = 5, Vanilla NMF")
nmf_k8 <- nmf_k8_p1 + ggtitle("K = 8, Vanilla NMF")

plot_grid(plotlist = c(list(nmf_k4), KimPark_k4_plots,
                       list(nmf_k5), KimPark_k5_plots,
                       list(nmf_k8), KimPark_k8_plots),
          nrow = 5, ncol = 3, byrow = FALSE)
ggsave("./output/plots/ex_KimPark.pdf", width = 7, height = 7)

plot_grid(plotlist = c(list(nmf_k4), Hoyer_k4_plots,
                       list(nmf_k5), Hoyer_k5_plots,
                       list(nmf_k8), Hoyer_k8_plots),
          nrow = 5, ncol = 3, byrow = FALSE)
ggsave("./output/plots/ex_Hoyer.pdf", width = 7, height = 7)


rcppml_k4_plots <- list()
rcppml_k4_t <- numeric()
for (L1pen in c(.01, 0.1, 0.3, 0.6)) {
  next_res <- run_RcppML_sparse_nmf(Y, k = 4, L1pen = L1pen, seeds = 1:10)
  next_p <- plot_nmf(next_res$fit, L, pops, paste0("K = 4, L1 penalty = ", L1pen), gap = 50)
  rcppml_k4_plots <- c(rcppml_k4_plots, list(next_p))
  rcppml_k4_t <- c(rcppml_k4_t, next_res$t[3])
}

rcppml_k5_plots <- list()
rcppml_k5_t <- numeric()
for (L1pen in c(.01, 0.1, 0.3, 0.6)) {
  next_res <- run_RcppML_sparse_nmf(Y, k = 5, L1pen = L1pen, seeds = 1:10)
  next_p <- plot_nmf(next_res$fit, L, pops, paste0("K = 5, L1 penalty = ", L1pen), gap = 50)
  rcppml_k5_plots <- c(rcppml_k5_plots, list(next_p))
  rcppml_k5_t <- c(rcppml_k5_t, next_res$t[3])
}

rcppml_k8_plots <- list()
rcppml_k8_t <- numeric()
for (L1pen in c(.01, 0.1, 0.3, 0.6)) {
  next_res <- run_RcppML_sparse_nmf(Y, k = 8, L1pen = L1pen, seeds = 1:10)
  next_p <- plot_nmf(next_res$fit, L, pops, paste0("K = 8, L1 penalty = ", L1pen), gap = 50)
  rcppml_k8_plots <- c(rcppml_k8_plots, list(next_p))
  rcppml_k8_t <- c(rcppml_k8_t, next_res$t[3])
}


plot_grid(plotlist = c(list(nmf_k4), rcppml_k4_plots,
                       list(nmf_k5), rcppml_k5_plots,
                       list(nmf_k8), rcppml_k8_plots),
          nrow = 5, ncol = 3, byrow = FALSE)
ggsave("output/plots/ex_rcppml.pdf", width = 7, height = 7)


ebnmf_res1 <- run_ebnmf_from_nmf(Y, nmf_k8_res$fit)
ebnmf_p1 <- plot_fl(ebnmf_res1$fit, L, pops, "EBNMF, NMF init (K = 8)", gap = 50)

ebnmf_res2 <- run_greedy_backfit(Y, Kmax = 10, alternating = FALSE)
ebnmf_p2 <- plot_fl(ebnmf_res2$fit, L, pops, "EBNMF, greedy/backfit", gap = 50)

plot_grid(ebnmf_p1, ebnmf_p2, nrow = 1)
ggsave("output/plots/ex_ebnmf.pdf", width = 7, height = 2.5)


t_df <- tibble(method = "KimPark", t = c(KimPark_k4_t, KimPark_k5_t, KimPark_k8_t)) |>
  bind_rows(tibble(method = "RcppML", t = c(nmf_k4_res$t[3], nmf_k5_res$t[3], nmf_k8_res$t[3], rcppml_k4_t, rcppml_k5_t, rcppml_k8_t))) |>
  bind_rows(tibble(method = "Hoyer", t = c(Hoyer_k4_t, Hoyer_k5_t, Hoyer_k8_t))) |>
  bind_rows(tibble(method = "NMF-EBNMF", t = ebnmf_res1$t[3])) |>
  bind_rows(tibble(method = "EBNMF-gb", t = ebnmf_res2$t[3])) |>
  mutate(method = factor(method, levels = c(
    "RcppML", "KimPark", "Hoyer", "NMF-EBNMF", "EBNMF-gb"
  )))
ggplot(t_df, aes(x = method, y = t)) +
  geom_boxplot() +
  scale_y_log10() +
  labs(x = "", y = "Elapsed time") +
  theme_minimal()
ggsave("./output/plots/ex_timing.pdf", width = 7, height = 2.5)


L1pens <- c(0, .001) #, .01, .05, seq(0.1, 0.8, by = 0.1))
k4_cv <- run_RcppML_cv(Y, k = 4, L1pens = L1pens, nfolds = 100, ntrials = 2, verbose = TRUE)
k5_cv <- run_RcppML_cv(Y, k = 5, L1pens = L1pens, nfolds = 100, ntrials = 100, verbose = TRUE)
k8_cv <- run_RcppML_cv(Y, k = 8, L1pens = L1pens, nfolds = 100, ntrials = 100, verbose = TRUE)
cv_res <- k4_cv |> mutate(k = "4") |>
  bind_rows(k5_cv |> mutate(k = "5")) #|>
  #bind_rows(k8_cv |> mutate(k = "8"))
saveRDS(cv_res, "./output/ex_cv.rds")

ggplot(cv_res, aes(x = L1pen, y = RMSE, linetype = k, color = k)) +
  geom_line() +
  scale_color_brewer(palette = "Purples") +
  theme_minimal() +
  labs(x = "L1 penalty", y = "RMSE", linetype = "K")
ggsave("output/plots/cv_rcppml.pdf", width = 5, height = 3)
