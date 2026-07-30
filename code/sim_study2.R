source("./code/sim_functions.R")

matlab_bin <- "/Applications/MATLAB_R2025a.app/bin/matlab"
options(matlab = paste(
  shQuote(matlab_bin),
  "-nodesktop -nosplash -nodisplay"
))

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
        bind_rows(run_all_methods(
          sim_dat, which_dat, Kmax, next_seed, colwise_noise, varied_n,
          ntrials = 10, verbose = verbose
        ))
    }
  }

  return(all_res)
}

setwd("./matlab/")
Matlab$startServer()
matlab <- Matlab()
open(matlab)

setVerbose(matlab, threshold = Inf)
# To suppress the annoying MATLAB chatter, edit inst/externals/MatlabServer.m
#   in the R.Matlab package.

sims_K3 <- run_sims(matlab, "Y", Kmax = 3, verbose = TRUE)
saveRDS(sims_K3, "../output/ss2_K3.rds")
sims_K4 <- run_sims(matlab, "Y", Kmax = 4, verbose = TRUE)
saveRDS(sims_K4, "../output/ss2_K4.rds")
sims_K6 <- run_sims(matlab, "Y", Kmax = 6, verbose = TRUE)
saveRDS(sims_K6, "../output/ss2_K6.rds")
sims_K8 <- run_sims(matlab, "Y", Kmax = 8, verbose = TRUE)
saveRDS(sims_K8, "../output/ss2_K8.rds")

sims_K3_colwise <- run_sims(matlab, "Y", Kmax = 3, verbose = TRUE, colwise_noise = TRUE)
saveRDS(sims_K3_colwise, "../output/ss2_K3_colwise.rds")
sims_K4_colwise <- run_sims(matlab, "Y", Kmax = 4, verbose = TRUE, colwise_noise = TRUE)
saveRDS(sims_K4_colwise, "../output/ss2_K4_colwise.rds")
sims_K6_colwise <- run_sims(matlab, "Y", Kmax = 6, verbose = TRUE, colwise_noise = TRUE)
saveRDS(sims_K6_colwise, "../output/ss2_K6_colwise.rds")
sims_K8_colwise <- run_sims(matlab, "Y", Kmax = 8, verbose = TRUE, colwise_noise = TRUE)
saveRDS(sims_K8_colwise, "../output/ss2_K8_colwise.rds")

close(matlab)
setwd("../")
