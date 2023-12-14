swimmer <- R.matlab::readMat("data/Y.mat")
dat <- apply(swimmer$Y, 3, as.vector)

dat <- dat - 1
rows_to_keep <- which(rowSums(dat) > 0)
dat <- dat[rows_to_keep, ]

t0 <- Sys.time()
fl <- flashier::flash(dat, ebnm_fn = ebnm_point_exponential, backfit = TRUE)
tf <- Sys.time() - t0

library(tidyverse)
fitted_L <- matrix(0, nrow = 32^2, ncol = fl$n_factors)
fitted_L[rows_to_keep, ] <- ldf(fl, type = "m")$L

tib <- as_tibble(fitted_L) |>
  mutate(
    row = rep(1:32, each = 32),
    col = rep(1:32, times = 32)
  ) |>
  pivot_longer(
    cols = -c(row, col),
    names_to = "k",
    values_to = "loading",
    names_prefix = "V",
    names_transform = as.numeric
  )
ggplot(tib, aes(x = row, y = col, fill = loading)) +
  geom_tile() +
  scale_y_reverse() +
  scale_fill_gradient(low = "white", high = "black") +
  facet_wrap(~k) +
  guides(fill = "none")
