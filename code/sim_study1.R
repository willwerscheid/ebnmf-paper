source("./code/sim_functions.R")
options(matlab = "/Applications/MATLAB_R2025a.app/bin/matlab")


### Simulations -----

sim_data <- function(ns, p, colwise_noise, component_scales = rep(4, length(ns))) {
  L <- matrix(0, nrow = sum(ns), ncol = length(ns))
  pops <- rep(LETTERS[1:length(ns)], times = ns)
  L[pops == "A", 1] <- 1
  L[pops == "B", 2] <- 1
  L[pops == "C", 3] <- 1

  L <- vary_library_sizes(L)
  F <- sim_F(p, ncol(L))
  F <- t(t(F) * component_scales)
  X <- sim_X(L, F, colwise_noise = colwise_noise)
  return_sim_data(X, L, F, pops)
}

run_sims <- function(matlab, which_dat, Kmax, nreps = 5, verbose = FALSE,
                     colwise_noise = FALSE) {
  all_res <- tibble()
  next_seed <- 0
  for (varied_n in c(2, 3, 5, 10, 15, 20, 25, 30)) {
    for (reps in 1:nreps) {
      if (verbose) cat("RARE N:", varied_n, " SEED:", next_seed + 1, "\n")

      ns <- c(1000, 250, varied_n)
      p <- 2000

      next_seed <- next_seed + 1
      set.seed(next_seed)
      sim_dat <- sim_data(ns, p, colwise_noise)

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
saveRDS(sims_K3, "../output/ss1_K3.rds")
sims_K4 <- run_sims(matlab, "Y", Kmax = 4, verbose = TRUE)
saveRDS(sims_K4, "../output/ss1_K4.rds")
sims_K6 <- run_sims(matlab, "Y", Kmax = 6, verbose = TRUE)
saveRDS(sims_K6, "../output/ss1_K6.rds")

# sims_K3 <- run_sims(matlab, "Ylibnorm", Kmax = 3, verbose = TRUE)
# saveRDS(sims_K3, "../output/ss1_K3_libnorm.rds")
# sims_K4 <- run_sims(matlab, "Ylibnorm", Kmax = 4, verbose = TRUE)
# saveRDS(sims_K4, "../output/ss1_K4_libnorm.rds")
# sims_K6 <- run_sims(matlab, "Ylibnorm", Kmax = 6, verbose = TRUE)
# saveRDS(sims_K6, "../output/ss1_K6_libnorm.rds")
#
# sims_K3 <- run_sims(matlab, "Y", Kmax = 3, verbose = TRUE, colwise_noise = TRUE)
# saveRDS(sims_K3, "../output/ss1_K3_colwise.rds")
# sims_K4 <- run_sims(matlab, "Y", Kmax = 4, verbose = TRUE, colwise_noise = TRUE)
# saveRDS(sims_K4, "../output/ss1_K4_colwise.rds")
# sims_K6 <- run_sims(matlab, "Y", Kmax = 6, verbose = TRUE, colwise_noise = TRUE)
# saveRDS(sims_K6, "../output/ss1_K6_colwise.rds")

close(matlab)
setwd("../")
