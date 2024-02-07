library(imager)
library(flashier)
library(tidyverse)
library(tictoc)

img_files <- list.files(
  path = "./data/faces/train/face/",
  pattern = ".*pgm",
  recursive = TRUE,
  full.names = TRUE
)

all_imgs <- lapply(img_files, pixmap::read.pnm)
all_imgs <- lapply(all_imgs, pixmap::getChannels)
all_imgs <- lapply(all_imgs, as.vector)

mat <- matrix(unlist(all_imgs), nrow = length(all_imgs[[1]]))
saveRDS(mat, "./data/cbcl.rds")

tic()
min_sd <- min(apply(tmp, 1, sd)) / 10
old_nfactors <- 0
fl <- flash_init(mat, var_type = 1, S = min_sd) |>
  flash_set_verbose(1) |>
  flash_greedy(Kmax = 40, ebnm_fn = ebnm_point_exponential)

while(fl$n_factors > old_nfactors) {
  old_nfactors <- fl$n_factors
  fl <- flash_backfit(fl, maxiter = 2)
  fl <- fl |>
    flash_greedy(Kmax = 40 - fl$n_factors, ebnm_fn = ebnm_point_exponential)
}

flb <- flash_backfit(fl, maxiter = 200, verbose = 3)
toc()

fitted_L <- ldf(fl, type = "m")$L

n <- 19
tib <- as_tibble(fitted_L) |>
  mutate(
    col = rep(1:n, times = n),
    row = rep(1:n, each = n)
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
  scale_fill_gradient(low = "black", high = "white") +
  facet_wrap(~k) +
  guides(fill = "none")
