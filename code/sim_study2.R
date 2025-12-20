source("./code/sim_functions.R")
options(matlab = "/Applications/MATLAB_R2025a.app/bin/matlab")


### Simulations -----

sim_data <- function(ns, p, colwise_noise, component_scales = rep(5, length(ns))) {
  L <- matrix(0, nrow = sum(ns), ncol = 3)
  pops <- rep(LETTERS[1:length(ns)], times = ns)
  L[pops == "A", 1] <- 1
  L[pops == "B", 2] <- 1
  L[pops == "C", 3] <- 1

  pi1 <- seq(0, 1, length.out = ns[4])
  pi2 <- seq(0, 1, length.out = ns[5])
  pi3 <- seq(0, 1, length.out = ns[6])
  L[pops == "D", 1:2] <- c(pi1, 1 - pi1)
  L[pops == "E", c(1, 3)] <- c(pi2, 1 - pi2)
  L[pops == "F", 2:3] <- c(pi3, 1 - pi3)

  L <- vary_library_sizes(L)
  F <- sim_F(p, ncol(L))
  F <- t(t(F) * component_scales)
  X <- sim_X(L, F, colwise_noise)
  return_sim_data(X, L, F, pops)
}

run_sims <- function(matlab, which_dat, Kmax, nreps = 5, verbose = FALSE,
                     colwise_noise = FALSE) {
  all_res <- tibble()
  next_seed <- 0
  for (varied_n in c(0, 1, 5, 50, 250, 450, 495, 499, 500)) {
    for (reps in 1:nreps) {
      if (verbose) cat("ANCHOR N:", varied_n, " SEED:", next_seed + 1, "\n")

      ns <- c(rep(varied_n, 3), rep(500 - varied_n, 3))
      p <- 2000

      next_seed <- next_seed + 1
      set.seed(next_seed)
      sim_dat <- sim_data(ns, p, colwise_noise = colwise_noise)

      all_res <- all_res |>
        bind_rows(run_all_methods(sim_dat, which_dat, Kmax, next_seed, colwise_noise, varied_n))
    }
  }

  return(all_res)
}

setwd("./matlab/")
Matlab$startServer(matlab = "/Applications/MATLAB_R2025a.app/bin/matlab")
matlab <- Matlab()
open(matlab)
setVerbose(matlab, threshold = 1000000)

sims_K3 <- run_sims(matlab, "Y", Kmax = 3, verbose = TRUE)
saveRDS(sims_K3, "../output/ss2_K3.rds")
sims_K4 <- run_sims(matlab, "Y", Kmax = 4, verbose = TRUE)
saveRDS(sims_K4, "../output/ss2_K4.rds")
sims_K6 <- run_sims(matlab, "Y", Kmax = 6, verbose = TRUE)
saveRDS(sims_K6, "../output/ss2_K6.rds")
sims_K8 <- run_sims(matlab, "Y", Kmax = 8, verbose = TRUE)
saveRDS(sims_K8, "../output/ss2_K8.rds")

sims_K3 <- run_sims(matlab, "Ylibnorm", Kmax = 3, verbose = TRUE)
saveRDS(sims_K3, "../output/ss2_K3_libnorm.rds")
sims_K4 <- run_sims(matlab, "Ylibnorm", Kmax = 4, verbose = TRUE)
saveRDS(sims_K4, "../output/ss2_K4_libnorm.rds")
sims_K6 <- run_sims(matlab, "Ylibnorm", Kmax = 6, verbose = TRUE)
saveRDS(sims_K6, "../output/ss2_K6_libnorm.rds")
sims_K8 <- run_sims(matlab, "Ylibnorm", Kmax = 8, verbose = TRUE)
saveRDS(sims_K8, "../output/ss2_K8_libnorm.rds")

# sims_K3 <- run_sims(matlab, "Y", Kmax = 3, verbose = TRUE, colwise_noise = TRUE)
# saveRDS(sims_K3, "../output/ss2_K3_colwise.rds")
# sims_K4 <- run_sims(matlab, "Y", Kmax = 4, verbose = TRUE, colwise_noise = TRUE)
# saveRDS(sims_K4, "../output/ss2_K4_colwise.rds")
# sims_K6 <- run_sims(matlab, "Y", Kmax = 6, verbose = TRUE, colwise_noise = TRUE)
# saveRDS(sims_K6, "../output/ss2_K6_colwise.rds")
# sims_K8 <- run_sims(matlab, "Y", Kmax = 8, verbose = TRUE, colwise_noise = TRUE)
# saveRDS(sims_K8, "../output/ss2_K8_colwise.rds")

close(matlab)
setwd("../")


### Plots -----

plot_df <- sims_K6 |> mutate(K = "K = 3") |>
  #bind_rows(sims_K4 |> mutate(K = "K = 4")) |>
  #bind_rows(sims_K6 |> mutate(K = "K = 6")) |>
  mutate(varied_n = shape) |>
  group_by(K, method, submethod, varied_n, seed) |>
  summarize(metric_val = max(metric_val * (metric_type == "FFcosine3"))) |>
  mutate(metric_val = coalesce(metric_val, 0)) |>
  mutate(varied_n = factor(varied_n)) |>
  mutate(method = fct_rev(method), submethod = fct_rev(submethod)) |>
  group_by(K, method, submethod, varied_n) |>
  summarize(metric_val = mean(metric_val)) |>
  ungroup()

all_plots <- list()
methods <- c("NMF", "EBNMF", "KimPark", "RcppML", "Hoyer")
margins <- c(80, 1, 74, 82, 82)
for (i in 1:length(methods)) {
  next_df <- plot_df |>
    filter(method == methods[[i]]) |>
    mutate(submethod = fct_drop(submethod))
  all_plots[[i]] <- ggplot(next_df, aes(x = varied_n, y = submethod, fill = metric_val)) +
    geom_tile() +
    scale_fill_gradient(low = "red4", high = "seagreen4", limits = c(0.95, 1)) +
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
p1 <- plot_grid(plotlist = all_plots, ncol = 1, rel_heights = 3 + c(3, 8, 9, 8, 9))
p2 <- plot_grid(p1, legend, nrow = 1, rel_widths = c(9, 1))
p3 <- add_sub(p2, label = "Size of \"pure\" populations",
              hjust = 0.4, vjust = -0.8, size = 10)
plot(p3)
#save_plot("./output/plots/ss2_FF1_allK.pdf", p3, base_width = 7, base_height = 5)
