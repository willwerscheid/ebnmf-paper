library(tidyverse)
library(cowplot)
library(ggrepel)
library(ggpubr)
library(reactable)
library(htmltools)
library("htmlwidgets")
library("webshot2")

glasbey <- fastTopics:::glasbey()[1:200]

scale_FF <- function(FF, LL, D = 1) {
  FF_norms <- apply(FF, 2, max)
  FF <- t(t(FF) / FF_norms)
  return(FF)
}
scale_LL <- function(FF, LL, D = 1) {
  FF_norms <- apply(FF, 2, max)
  LL <- t(t(LL) * FF_norms)
  return(LL)
}

ebnmf30 <- readRDS("output/montoro-ebmf.rds")
ebnmf30_FF <- scale_FF(ebnmf30$fit$F_pm, ebnmf30$fit$L_pm)

ebnmf40 <- readRDS("output/montoro-ebmf-40.rds")
ebnmf40_FF <- scale_FF(ebnmf40$fit$F_pm, ebnmf40$fit$L_pm)

ebnmf <- ebnmf30
ebnmf_FF <- ebnmf30_FF

cell_types <- c(
  "Ciliated",
  "Proliferating",
  "Basal",
  "Club",
  "Club (hillock-associated)",
  "Goblet.progenitor",
  "Goblet.1",
  "Goblet.2",
  "Tuft.progenitor",
  "Tuft.1",
  "Tuft.2",
  "Neuroendocrine",
  "Ionocyte"
)

ct_abbr <- cell_types
ct_abbr <- str_replace(ct_abbr, "\\.", "\\-")
ct_abbr[2] <- "Prolif."
ct_abbr[5] <- "Hillock"
ct_abbr[12] <- "NEC"
ct_abbr <- str_replace(ct_abbr, "progenitor", "Pr.")

external_info <- tibble(
  CellIdx = 1:nrow(ebnmf30_FF),
  CellType = sapply(strsplit(rownames(ebnmf30_FF), "_"), `[[`, 5),
  TimePoint = sapply(strsplit(rownames(ebnmf30_FF), "_"), `[[`, 2),
  Replicate = sapply(strsplit(rownames(ebnmf30_FF), "_"), `[[`, 3)
) |>
  left_join(tibble(CellType = cell_types, CtAbbr = ct_abbr), by = "CellType") |>
  mutate(
    CellType = factor(CellType, levels = cell_types),
    CtAbbr = factor(CtAbbr, levels = ct_abbr),
    TimePoint = factor(TimePoint),
    Replicate = factor(Replicate),
  ) |>
  mutate(
    TPxRep = fct_cross(TimePoint, Replicate),
  )


# Make data frames. -----

make_tib <- function(FF, sample_idx, info_type = "CtAbbr") {
  info <- external_info[[info_type]]

  FF <- FF[sample_idx, ]
  info <- info[sample_idx]
  colnames(FF) <- paste0("k", 1:ncol(FF))

  tib <- as_tibble(FF) |>
    mutate(
      CellIdx = row_number(),
      Info = info
    ) |>
    pivot_longer(
      -c(CellIdx, Info),
      names_to = "Component",
      values_to = "Loading",
      values_drop_na = TRUE
    ) |>
    group_by(CellIdx, Info, Component) |>
    summarize(Loading = sum(Loading), .groups = "drop") |>
    mutate(
      Component = factor(Component, levels = colnames(FF))
    ) |>
    arrange(CellIdx, Info, Component)

  return(tib)
}


# Order loadings. -----

do_order <- function(tib, seed = 1) {
  set.seed(seed)

  idx <- numeric(0)
  for (grp in levels(tib$Info)) {
    next_grp <- tib |>
      filter(Info == grp) |>
      pivot_wider(names_from = Component, values_from = Loading)
    Lmat <- next_grp |>
      select(-(1:2)) |>
      as.matrix()

    Fmat <- matrix(1, nrow(Lmat), ncol(Lmat))
    fit <- list(F = Fmat, L = Lmat)
    class(fit) <- "poisson_nmf_fit"

    tsne_res <- fastTopics::tsne_from_topics(
      fit, dims = 1, max_iter = 2500
    )

    idx <- c(idx, next_grp$CellIdx[order(tsne_res)])
  }

  return(idx)
}


# Stacked bar plots. -----

do_plot <- function(tib, loadings_order, colors, ptitle,
                    legend_ncol = 1) {
  tib <- tib |>
    pivot_wider(names_from = Component, values_from = Loading)
  group <- tib$Info

  Lmat <- tib |>
    select(-(1:2)) |>
    as.matrix()

  if (length(unique(group)) == 1) {
    p <- fastTopics::structure_plot(
      Lmat, topics = colnames(Lmat),
      loadings_order = rank(loadings_order),
      colors = colors,
      gap = 10,
      verbose = FALSE
    )
  } else {
    p <- fastTopics::structure_plot(
      Lmat, topics = colnames(Lmat),
      grouping = group,
      loadings_order = rank(loadings_order),
      colors = colors,
      gap = 10,
      verbose = FALSE
    )
  }

  p <- p + labs(
    y = "Membership",
    title = ptitle,
    fill = "Component"
  ) +
    guides(
      fill = guide_legend(ncol = legend_ncol)
    ) +
    theme(legend.position = "right")

  return(p)
}


# Top genes. -----

ebnmf_LL <- scale_LL(ebnmf$fit$F_pm, ebnmf$fit$L_pm)
colnames(ebnmf_LL) <- paste0("k", 1:ncol(ebnmf_LL))
ebnmf_LL_tib <- as_tibble(ebnmf_LL) |>
  mutate(Gene = rownames(ebnmf_LL)) |>
  pivot_longer(-Gene, names_to = "Component", values_to = "Loading") |>
  mutate(Component = factor(Component, levels = colnames(ebnmf_LL))) |>
  group_by(Gene) |>
  mutate(Uniqueness = Loading / sort(Loading, decreasing = TRUE)[2])

top_unique_genes <- ebnmf_LL_tib |>
  filter(Loading > 0.5, Uniqueness > 1) |>
  mutate(Gene = str_remove(Gene, "_.*")) |>
  group_by(Component) |>
  slice_max(Uniqueness, n = 15) |>
  summarize(TopGenes = paste(Gene, collapse = ", "))

top_unique_genes <- top_unique_genes |>
  mutate(ComponentName = paste0(Component, " (", TopGenes, ")"))


# Main EBNMF figure. -----

# set.seed(1)
# all_samp <- sample(1:66265, 6000)
# ebnmf_all <- make_tib(ebnmf_FF, all_samp) |>
#   mutate(Info = factor("Any"))
# ebnmf_all_idx <- do_order(ebnmf_all)
# ebnmf_all_p <- do_plot(ebnmf_all, ebnmf_all_idx, glasbey[2:41],
#                        "All cells", legend_ncol = 2)

set.seed(1)
ct_samp <- external_info |>
  slice_sample(n = 200, by = CellType, replace = TRUE) |>
  pull(CellIdx)
ebnmf_ct <- make_tib(ebnmf_FF, ct_samp)
ebnmf_ct_idx <- do_order(ebnmf_ct)

plot_df <- ebnmf_ct |>
  left_join(top_unique_genes, by = "Component") |>
  mutate(Component = factor(
    ComponentName, levels = top_unique_genes$ComponentName)
  ) |>
  select(-TopGenes, -ComponentName)

ebnmf_ct_p <- do_plot(plot_df, ebnmf_ct_idx, glasbey[2:41],
                      "EBNMF results", legend_ncol = 1) +
  labs(fill = "Component (Top unique genes)") +
  theme(legend.position = "bottom",
        legend.title.position = "top",
        legend.title = element_text(face = "bold"))

# ebnmf_legend <- as_ggplot(get_legend(ebnmf_ct_p))
#
# p <- plot_grid(
#   ebnmf_ct_p + theme(legend.position = "none"),
#   ebnmf_legend,
#   nrow = 2, ncol = 1,
#   rel_heights = c(4, 6)
# )

save_plot("./output/plots/montoro_ebnmf.pdf",
          ebnmf_ct_p, base_height = 10, base_width = 8)




nmf <- readRDS("output/montoro-nmf.rds")
nmf_FF <- scale_FF(t(nmf$fit$h), nmf$fit$w, nmf$fit$d)

nmf40 <- readRDS("output/montoro-nmf-40.rds")
nmf40_FF <- scale_FF(t(nmf40$fit$h), nmf40$fit$w, nmf40$fit$d)

tm <- readRDS("output/montoro-topics.rds")
tm_fit <- tm$fit
tm_fit$F <- tm$fit$L
tm_fit$L <- tm$fit$F
tm_fit$Fn <- tm$fit$Ln
tm_fit$Ln <- tm$fit$Fn
tm_fit$Fy <- tm$fit$Ly
tm_fit$Ly <- tm$fit$Fy
tm_fit <- fastTopics::poisson2multinom(tm_fit)
tm_FF <- tm_fit$L

tm40 <- readRDS("output/montoro-topics-40.rds")
tm_fit <- tm40$fit
tm_fit$F <- tm40$fit$L
tm_fit$L <- tm40$fit$F
tm_fit$Fn <- tm40$fit$Ln
tm_fit$Ln <- tm40$fit$Fn
tm_fit$Fy <- tm40$fit$Ly
tm_fit$Ly <- tm40$fit$Fy
tm_fit <- fastTopics::poisson2multinom(tm_fit)
tm40_FF <- tm_fit$L




# Select relevant components. -----

select_k <- function(FF, sample_idx, which_types, unif_sparse_k = NULL) {
  y <- external_info |>
    # filter(CellIdx %in% sample_idx) |>
    mutate(Y = 1L * str_detect(CtAbbr, which_types)) |>
    pull(Y)
  coefs <- numeric(30)
  for (k in 1:30) {
    if (k %in% unif_sparse_k) {
      coefs[k] <- -Inf
    } else {
      lmod <- glm(y ~ FF[, k], family = "binomial")
      #lmod <- glm(y ~ FF[sample_idx, k], family = "binomial")
      if (summary(lmod)$coef[2, 4] > .1) {
        coefs[k] <- 0
      } else {
        coefs[k] <- coef(lmod)[2]
      }
    }
  }
  which_k <- order(coefs, decreasing = TRUE)[1:sum(coefs > 0)]
  return(rev(which_k))
}

select_k_tp <- function(FF, sample_idx, which_tp, unif_sparse_k = NULL) {
  y <- external_info |>
    filter(CellIdx %in% sample_idx) |>
    mutate(Y = 1L * str_detect(TimePoint, which_tp)) |>
    pull(Y)
  coefs <- numeric(30)
  for (k in 1:30) {
    if (k %in% unif_sparse_k) {
      coefs[k] <- -Inf
    } else {
      lmod <- glm(y ~ FF[sample_idx, k], family = "binomial")
      if (summary(lmod)$coef[2, 4] > .1) {
        coefs[k] <- 0
      } else {
        coefs[k] <- coef(lmod)[2]
      }
    }
  }
  which_k <- order(coefs, decreasing = TRUE)[1:sum(coefs > 0)]
  return(rev(which_k))
}








# Structure plots. -----

do_plot <- function(tib, loadings_order, colors, ptitle, ntiles = 400) {
  tib <- tib |>
    filter(Component != "Other") |>
    mutate(Component = paste0("k", Component))

  # if (!is.null(tib$Label)) {
  #   tib <- tib |>
  #     mutate(Component = Label) |>
  #     select(-Label)
  # }

  struct_df <- tib |>
    pivot_wider(names_from = Component, values_from = Loading)
  Lmat <- struct_df |>
    select(-(1:2)) |>
    as.matrix()
  group <- struct_df$Info

  if ("Other" %in% tib$Component) {
    colors <- c("grey90", colors)
  }

  p <- fastTopics::structure_plot(
    Lmat, topics = colnames(Lmat),
    grouping = group,
    loadings_order = loadings_order,
    colors = colors,
    gap = 10,
    verbose = FALSE
  # ) + theme(
  #   axis.text.y = element_blank()
  ) + labs(
    y = "Membership",
    title = ptitle,
    fill = "Component"
  ) +
    guides(
      fill = guide_legend(ncol = 1)
    ) +
    theme(legend.position = "right", legend.title = element_blank())
    # guides(fill = "none")

  return(p)
}


# All cells (ungrouped). -----




# All types. -----

set.seed(1)
sample_idx <- external_info |>
  group_by(CellType) |>
  mutate(rand_n = sample(1:n())) |>
  mutate(keep_row = (rand_n <= 100)) |>
  filter(keep_row) |>
  pull(CellIdx)

ebnmf_tib <- make_tib(ebnmf_FF, sample_idx, 1:40)

combined_tib <- ebnmf_tib |> mutate(Method = "EBNMF")
set.seed(1)
idx <- numeric(0)
for (grp in levels(combined_tib$Info)) {
  idx <- c(idx, do_order(combined_tib, grp))
}

ebnmf_p <- do_plot(ebnmf_tib, idx, glasbey[2:41], "EBNMF")


# Tuft/NEC/Ionocytes. -----

cell_types <- "Tuf|NEC|Ion"

set.seed(1)
sample_idx <- external_info |>
  group_by(CellType) |>
  mutate(rand_n = sample(1:n())) |>
  mutate(keep_row = (rand_n <= 50 |
                       (str_detect(CtAbbr, cell_types) & rand_n <= 100))) |>
  filter(keep_row) |>
  pull(CellIdx)

ebnmf_k <- select_k(ebnmf_FF, sample_idx, cell_types)
nmf_k <- select_k(nmf_FF, sample_idx, cell_types) # 1, 26, 29
tm_k <- select_k(tm_FF, sample_idx, cell_types) # 7, 24, 2

ebnmf_tib <- make_tib(ebnmf_FF, sample_idx, ebnmf_k)
nmf_tib <- make_tib(nmf_FF, sample_idx, nmf_k)
tm_tib <- make_tib(tm_FF, sample_idx, tm_k)

combined_tib <- ebnmf_tib |> mutate(Method = "EBNMF") |>
  bind_rows(nmf_tib |> mutate(Method = "NMF")) |>
  bind_rows(tm_tib |> mutate(Method = "TM"))

set.seed(1)
idx <- numeric(0)
for (grp in levels(combined_tib$Info)) {
  idx <- c(idx, do_order(combined_tib, grp))
}

# ebnmf_tib <- ebnmf_tib |>
#   left_join(tibble(
#     Component = c("Other", "21", "30", "26", "19", "15", "12", "17"),
#     Label = c("X", "Ubiq.", "Sprs-a", "Sprs-b", "Ion.", "Tuft-1", "NEC", "Tuft-2")
#   ))
# nmf_tib <- nmf_tib |>
#   left_join(tibble(
#     Component = c("Other", "1", "26", "29"),
#     Label = c("X", "Ubiq.", "NEC/Ion.", "Tuft")
#   ))
# tm_tib <- tm_tib |>
#   left_join(tibble(
#     Component = c("Other", "7", "24", "2"),
#     Label = c("X", "Ubiq-a", "Ubiq-b", "G2/T/N/I")
#   ))

ebnmf_p <- do_plot(ebnmf_tib, idx, glasbey[2:8], "EBNMF") +
  theme(axis.text.x = element_blank())
nmf_p <- do_plot(nmf_tib, idx, glasbey[34:36], "NMF") +
  theme(axis.text.x = element_blank())
tm_p <- do_plot(tm_tib, idx, glasbey[61:64], "Topic Model")
tni_p <- plot_grid(ebnmf_p, nmf_p, tm_p, nrow = 3, ncol = 1,
                   rel_heights = c(3, 3, 3.9))
save_plot("./output/plots/montoro_tni.pdf", tni_p, base_height = 6.25, base_width = 8)


# Goblet. -----

cell_types <- "Gob"

set.seed(1)
sample_idx <- external_info |>
  group_by(CellType) |>
  mutate(rand_n = sample(1:n())) |>
  mutate(keep_row = (rand_n <= 30 |
                       (str_detect(CtAbbr, cell_types) & rand_n <= 100))) |>
  filter(keep_row) |>
  pull(CellIdx)

ebnmf_k <- select_k(ebnmf_FF, sample_idx, cell_types) # 3, 13, 29, 28, 24
nmf_k <- select_k(nmf_FF, sample_idx, cell_types) # 23, 26x, 1x, 19, 11, 25
tm_k <- select_k(tm_FF, sample_idx, cell_types) # 27, 25, 18, 16, 2, 13, 19, 7x, 26, 3, 17, 10

ebnmf_tib <- make_tib(ebnmf_FF, sample_idx, ebnmf_k)
nmf_tib <- make_tib(nmf_FF, sample_idx, nmf_k)
tm_tib <- make_tib(tm_FF, sample_idx, tm_k)

combined_tib <- ebnmf_tib |> mutate(Method = "EBNMF") |>
  bind_rows(nmf_tib |> mutate(Method = "NMF")) |>
  bind_rows(tm_tib |> mutate(Method = "TM"))

set.seed(1)
idx <- numeric(0)
for (grp in levels(combined_tib$Info)) {
  idx <- c(idx, do_order(combined_tib, grp))
}

# ebnmf_tib <- ebnmf_tib |>
#   left_join(tibble(
#     Component = c("Other", "3", "13", "29", "28", "24"),
#     Label = c("X", "Cl/Gob", "Gob/Cl", "Gob-2", "Gob-1sub", "Gob")
#   ))
# nmf_tib <- nmf_tib |>
#   left_join(tibble(
#     Component = c("Other", "23", "26", "1", "19", "11", "25"),
#     Label = c("X", "Cil.", "NEC/Ion.", "Ubiq.", "non-T/N/I", "Gob/Cl", "Gob")
#   ))
# tm_tib <- tm_tib |>
#   left_join(tibble(
#     Component = c("Other", "24", "2"),
#     Label = c("X", "Dense", "G2/T/N/I")
#   ))

ebnmf_p <- do_plot(ebnmf_tib, idx, glasbey[9:13], "EBNMF")
nmf_p <- do_plot(nmf_tib, idx, glasbey[40:46], "NMF")
tm_p <- do_plot(tm_tib, idx, glasbey[64:75], "Topic Model")
gob_p <- plot_grid(ebnmf_p, nmf_p, tm_p, nrow = 3, ncol = 1)
save_plot("./output/plots/montoro_gob.pdf", gob_p, base_height = 6.25, base_width = 8)


# Basal. -----

cell_types <- "Bas"

set.seed(1)
sample_idx <- external_info |>
  group_by(CellType) |>
  mutate(rand_n = sample(1:n())) |>
  mutate(keep_row = (rand_n <= 50 |
                       (str_detect(CtAbbr, cell_types) & rand_n <= 600))) |>
  filter(keep_row) |>
  pull(CellIdx)

ebnmf_k <- select_k(ebnmf_FF, sample_idx, cell_types) # 21x, 27, 7, 2, 6, 4, 9, 1
nmf_k <- select_k(nmf_FF, sample_idx, cell_types) # 12, 30, 24, 1, 13, 20
tm_k <- select_k(tm_FF, sample_idx, cell_types)

ebnmf_tib <- make_tib(ebnmf_FF, sample_idx, ebnmf_k)
nmf_tib <- make_tib(nmf_FF, sample_idx, nmf_k)
tm_tib <- make_tib(tm_FF, sample_idx, tm_k)

combined_tib <- ebnmf_tib |> mutate(Method = "EBNMF") |>
  bind_rows(nmf_tib |> mutate(Method = "NMF")) |>
  bind_rows(tm_tib |> mutate(Method = "TM"))

set.seed(1)
idx <- numeric(0)
for (grp in levels(combined_tib$Info)) {
  idx <- c(idx, do_order(combined_tib, grp))
}

ebnmf_p <- do_plot(ebnmf_tib, idx, glasbey[14:21], "EBNMF") +
  theme(axis.text.x = element_blank())
nmf_p <- do_plot(nmf_tib, idx, glasbey[47:62], "NMF") +
  theme(axis.text.x = element_blank())
tm_p <- do_plot(tm_tib, idx, glasbey[76:89], "Topic Model")
bas_p <- plot_grid(ebnmf_p, nmf_p, tm_p, nrow = 3, ncol = 1,
                   rel_heights = c(4.5, 4.5, 5.5))
save_plot("./output/plots/montoro_bas.pdf", bas_p, base_height = 7.5, base_width = 8)


# Ciliated. -----

cell_types <- "Cil"

set.seed(1)
sample_idx <- external_info |>
  group_by(CellType) |>
  mutate(rand_n = sample(1:n())) |>
  mutate(keep_row = (rand_n <= 50 |
                       (str_detect(CtAbbr, cell_types) & rand_n <= 600))) |>
  filter(keep_row) |>
  pull(CellIdx)

ebnmf_k <- select_k(ebnmf_FF, sample_idx, cell_types) # 21x, 27, 7, 2, 6, 4, 9, 1
nmf_k <- select_k(nmf_FF, sample_idx, cell_types) # 12, 30, 24, 1, 13, 20
tm_k <- select_k(tm_FF, sample_idx, cell_types)

ebnmf_tib <- make_tib(ebnmf_FF, sample_idx, ebnmf_k)
nmf_tib <- make_tib(nmf_FF, sample_idx, nmf_k)
tm_tib <- make_tib(tm_FF, sample_idx, tm_k)

combined_tib <- ebnmf_tib |> mutate(Method = "EBNMF") |>
  bind_rows(nmf_tib |> mutate(Method = "NMF")) |>
  bind_rows(tm_tib |> mutate(Method = "TM"))

set.seed(1)
idx <- numeric(0)
for (grp in levels(combined_tib$Info)) {
  idx <- c(idx, do_order(combined_tib, grp))
}

ebnmf_p <- do_plot(ebnmf_tib, idx, glasbey[2:10], "EBNMF") +
  theme(axis.text.x = element_blank())
nmf_p <- do_plot(nmf_tib, idx, glasbey[34:50], "NMF") +
  theme(axis.text.x = element_blank())
tm_p <- do_plot(tm_tib, idx, glasbey[61:80], "Topic Model")
cil_p <- plot_grid(ebnmf_p, nmf_p, tm_p, nrow = 3, ncol = 1,
                   rel_heights = c(4.5, 4.5, 5.5))
save_plot("./output/plots/montoro_cil.pdf", cil_p, base_height = 7.5, base_width = 8)


# Proliferating. -----

cell_types <- "Pro"

set.seed(1)
sample_idx <- external_info |>
  group_by(CellType) |>
  mutate(rand_n = sample(1:n())) |>
  mutate(keep_row = (rand_n <= 50 |
                       (str_detect(CtAbbr, cell_types) & rand_n <= 600))) |>
  filter(keep_row) |>
  pull(CellIdx)

ebnmf_k <- select_k(ebnmf_FF, sample_idx, cell_types) # 21x, 27, 7, 2, 6, 4, 9, 1
nmf_k <- select_k(nmf_FF, sample_idx, cell_types) # 12, 30, 24, 1, 13, 20
tm_k <- select_k(tm_FF, sample_idx, cell_types)

ebnmf_tib <- make_tib(ebnmf_FF, sample_idx, ebnmf_k)
nmf_tib <- make_tib(nmf_FF, sample_idx, nmf_k)
tm_tib <- make_tib(tm_FF, sample_idx, tm_k)

combined_tib <- ebnmf_tib |> mutate(Method = "EBNMF") |>
  bind_rows(nmf_tib |> mutate(Method = "NMF")) |>
  bind_rows(tm_tib |> mutate(Method = "TM"))

set.seed(1)
idx <- numeric(0)
for (grp in levels(combined_tib$Info)) {
  idx <- c(idx, do_order(combined_tib, grp))
}

ebnmf_p <- do_plot(ebnmf_tib, idx, glasbey[2:10], "EBNMF") +
  theme(axis.text.x = element_blank())
nmf_p <- do_plot(nmf_tib, idx, glasbey[34:50], "NMF") +
  theme(axis.text.x = element_blank())
tm_p <- do_plot(tm_tib, idx, glasbey[61:80], "Topic Model")
pro_p <- plot_grid(ebnmf_p, nmf_p, tm_p, nrow = 3, ncol = 1,
                   rel_heights = c(4.5, 4.5, 5.5))
save_plot("./output/plots/montoro_pro.pdf", pro_p, base_height = 7.5, base_width = 8)


# Club. -----

cell_types <- "Club"

set.seed(1)
sample_idx <- external_info |>
  group_by(CellType) |>
  mutate(rand_n = sample(1:n())) |>
  mutate(keep_row = (rand_n <= 50 |
                       (str_detect(CtAbbr, cell_types) & rand_n <= 600))) |>
  filter(keep_row) |>
  pull(CellIdx)

ebnmf_k <- select_k(ebnmf_FF, sample_idx, cell_types) # 21x, 27, 7, 2, 6, 4, 9, 1
nmf_k <- select_k(nmf_FF, sample_idx, cell_types) # 12, 30, 24, 1, 13, 20
tm_k <- select_k(tm_FF, sample_idx, cell_types)

ebnmf_tib <- make_tib(ebnmf_FF, sample_idx, ebnmf_k)
nmf_tib <- make_tib(nmf_FF, sample_idx, nmf_k)
tm_tib <- make_tib(tm_FF, sample_idx, tm_k)

combined_tib <- ebnmf_tib |> mutate(Method = "EBNMF") |>
  bind_rows(nmf_tib |> mutate(Method = "NMF")) |>
  bind_rows(tm_tib |> mutate(Method = "TM"))

set.seed(1)
idx <- numeric(0)
for (grp in levels(combined_tib$Info)) {
  idx <- c(idx, do_order(combined_tib, grp))
}

ebnmf_p <- do_plot(ebnmf_tib, idx, glasbey[2:20], "EBNMF") +
  theme(axis.text.x = element_blank())
nmf_p <- do_plot(nmf_tib, idx, glasbey[34:50], "NMF") +
  theme(axis.text.x = element_blank())
tm_p <- do_plot(tm_tib, idx, glasbey[61:80], "Topic Model")
club_p <- plot_grid(ebnmf_p, nmf_p, tm_p, nrow = 3, ncol = 1,
                   rel_heights = c(4.5, 4.5, 5.5))
save_plot("./output/plots/montoro_club.pdf", club_p, base_height = 7.5, base_width = 8)


# Hillock. -----

cell_types <- "Hil"

set.seed(1)
sample_idx <- external_info |>
  group_by(CellType) |>
  mutate(rand_n = sample(1:n())) |>
  mutate(keep_row = (rand_n <= 50 |
                       (str_detect(CtAbbr, cell_types) & rand_n <= 600))) |>
  filter(keep_row) |>
  pull(CellIdx)

ebnmf_k <- select_k(ebnmf_FF, sample_idx, cell_types) # 21x, 27, 7, 2, 6, 4, 9, 1
nmf_k <- select_k(nmf_FF, sample_idx, cell_types) # 12, 30, 24, 1, 13, 20
tm_k <- select_k(tm_FF, sample_idx, cell_types)

ebnmf_tib <- make_tib(ebnmf_FF, sample_idx, ebnmf_k)
nmf_tib <- make_tib(nmf_FF, sample_idx, nmf_k)
tm_tib <- make_tib(tm_FF, sample_idx, tm_k)

combined_tib <- ebnmf_tib |> mutate(Method = "EBNMF") |>
  bind_rows(nmf_tib |> mutate(Method = "NMF")) |>
  bind_rows(tm_tib |> mutate(Method = "TM"))

set.seed(1)
idx <- numeric(0)
for (grp in levels(combined_tib$Info)) {
  idx <- c(idx, do_order(combined_tib, grp))
}

ebnmf_p <- do_plot(ebnmf_tib, idx, glasbey[2:20], "EBNMF") +
  theme(axis.text.x = element_blank())
nmf_p <- do_plot(nmf_tib, idx, glasbey[34:50], "NMF") +
  theme(axis.text.x = element_blank())
tm_p <- do_plot(tm_tib, idx, glasbey[61:80], "Topic Model")
hil_p <- plot_grid(ebnmf_p, nmf_p, tm_p, nrow = 3, ncol = 1,
                    rel_heights = c(4.5, 4.5, 5.5))
save_plot("./output/plots/montoro_hil.pdf", hil_p, base_height = 7.5, base_width = 8)


# Time points. -----

set.seed(1)
sample_idx <- external_info |>
  group_by(TPxRep) |>
  slice_sample(n = 200) |>
  pull(CellIdx)

tp <- "Tp0"

ebnmf_k <- select_k_tp(ebnmf_FF, sample_idx, tp)
nmf_k <- select_k_tp(nmf_FF, sample_idx, tp)
tm_k <- select_k_tp(tm_FF, sample_idx, tp)

ebnmf_tib <- make_tib(ebnmf_FF, sample_idx, ebnmf_k, info_type = "TPxRep")
nmf_tib <- make_tib(nmf_FF, sample_idx, nmf_k, info_type = "TPxRep")
tm_tib <- make_tib(tm_FF, sample_idx, tm_k, info_type = "TPxRep")

combined_tib <- ebnmf_tib |> mutate(Method = "EBNMF") |>
  bind_rows(nmf_tib |> mutate(Method = "NMF")) |>
  bind_rows(tm_tib |> mutate(Method = "TM"))

set.seed(1)
idx <- numeric(0)
for (grp in levels(combined_tib$Info)) {
  idx <- c(idx, do_order(combined_tib, grp))
}

ebnmf_p <- do_plot(ebnmf_tib, idx, glasbey[99:106], "EBNMF")
nmf_p <- do_plot(nmf_tib, idx, glasbey[108:117], "NMF")
tm_p <- do_plot(tm_tib, idx, glasbey[118:129], "Topic Model")
tp0_p <- plot_grid(ebnmf_p, nmf_p, tm_p, nrow = 3, ncol = 1)


# 40-factor T/N/I. -----

cell_types <- "Tuf|NEC|Ion"

set.seed(1)
sample_idx <- external_info |>
  group_by(CellType) |>
  mutate(rand_n = sample(1:n())) |>
  mutate(keep_row = (rand_n <= 50 |
                       (str_detect(CtAbbr, cell_types) & rand_n <= 100))) |>
  filter(keep_row) |>
  pull(CellIdx)

ebnmf_k <- select_k(ebnmf_FF, sample_idx, cell_types)
ebnmf_k <- select_k(ebnmf_FF, sample_idx, cell_types)
nmf_k <- select_k(nmf_FF, sample_idx, cell_types)
nmf40_k <- select_k(nmf40_FF, sample_idx, cell_types)
tm_k <- select_k(tm_FF, sample_idx, cell_types)
tm40_k <- select_k(tm40_FF, sample_idx, cell_types)

ebnmf_tib <- make_tib(ebnmf_FF, sample_idx, ebnmf_k)
ebnmf_tib <- make_tib(ebnmf_FF, sample_idx, ebnmf_k)
nmf_tib <- make_tib(nmf_FF, sample_idx, nmf_k)
nmf40_tib <- make_tib(nmf40_FF, sample_idx, nmf40_k)
tm_tib <- make_tib(tm_FF, sample_idx, tm_k)
tm40_tib <- make_tib(tm40_FF, sample_idx, tm40_k)

combined_tib <- ebnmf_tib |> mutate(Method = "EBNMF") |>
  bind_rows(ebnmf_tib |> mutate(Method = "EBNMF-40")) |>
  bind_rows(nmf_tib |> mutate(Method = "NMF")) |>
  bind_rows(nmf40_tib |> mutate(Method = "NMF-40")) |>
  bind_rows(tm_tib |> mutate(Method = "TM")) |>
  bind_rows(tm40_tib |> mutate(Method = "TM-40"))

set.seed(1)
idx <- numeric(0)
for (grp in levels(combined_tib$Info)) {
  idx <- c(idx, do_order(combined_tib, grp))
}

ebnmf_p <- do_plot(ebnmf_tib, idx, glasbey[2:10], "EBNMF (30 components)") +
  theme(axis.text.x = element_blank())
ebnmf_p <- do_plot(ebnmf_tib, idx, glasbey[c(2, 5:10)], "EBNMF (40 components)") +
  theme(axis.text.x = element_blank())
nmf_p <- do_plot(nmf_tib, idx, glasbey[34:40], "NMF (30 components)") +
  theme(axis.text.x = element_blank())
nmf40_p <- do_plot(nmf40_tib, idx, glasbey[34:40], "NMF (40 components)") +
  theme(axis.text.x = element_blank())
tm_p <- do_plot(ebnmf_tib, idx, glasbey[61:70], "Topic Model (30 components)") +
  theme(axis.text.x = element_blank())
tm40_p <- do_plot(ebnmf_tib, idx, glasbey[61:70], "Topic Model (40 components)")
tni40_p <- plot_grid(ebnmf_p, ebnmf_p, nmf_p, nmf40_p, tm_p, tm40_p,
                     nrow = 6, ncol = 1,
                   rel_heights = c(3, 3, 3, 3, 3, 3.9))
save_plot("./output/plots/montoro_tni40.pdf", tni40_p, base_height = 12.25, base_width = 8)


# 40-factor basal. -----

cell_types <- "Club"

set.seed(1)
sample_idx <- external_info |>
  group_by(CellType) |>
  mutate(rand_n = sample(1:n())) |>
  mutate(keep_row = (rand_n <= 50 |
                       (str_detect(CtAbbr, cell_types) & rand_n <= 600))) |>
  filter(keep_row) |>
  pull(CellIdx)

ebnmf_k <- select_k(ebnmf_FF, sample_idx, cell_types)
ebnmf_k <- select_k(ebnmf_FF, sample_idx, cell_types)
nmf_k <- select_k(nmf_FF, sample_idx, cell_types)
nmf40_k <- select_k(nmf40_FF, sample_idx, cell_types)
tm_k <- select_k(tm_FF, sample_idx, cell_types)
tm40_k <- select_k(tm40_FF, sample_idx, cell_types)

ebnmf_tib <- make_tib(ebnmf_FF, sample_idx, ebnmf_k)
ebnmf_tib <- make_tib(ebnmf_FF, sample_idx, ebnmf_k)
nmf_tib <- make_tib(nmf_FF, sample_idx, nmf_k)
nmf40_tib <- make_tib(nmf40_FF, sample_idx, nmf40_k)
tm_tib <- make_tib(tm_FF, sample_idx, tm_k)
tm40_tib <- make_tib(tm40_FF, sample_idx, tm40_k)

combined_tib <- ebnmf_tib |> mutate(Method = "EBNMF") |>
  bind_rows(ebnmf_tib |> mutate(Method = "EBNMF-40")) |>
  bind_rows(nmf_tib |> mutate(Method = "NMF")) |>
  bind_rows(nmf40_tib |> mutate(Method = "NMF-40")) |>
  bind_rows(tm_tib |> mutate(Method = "TM")) |>
  bind_rows(tm40_tib |> mutate(Method = "TM-40"))

set.seed(1)
idx <- numeric(0)
for (grp in levels(combined_tib$Info)) {
  idx <- c(idx, do_order(combined_tib, grp))
}

ebnmf_p <- do_plot(ebnmf_tib, idx, glasbey[2:11], "EBNMF (30 components)") +
  theme(axis.text.x = element_blank()) +
  guides(
    fill = guide_legend(ncol = 2)
  )
ebnmf_p <- do_plot(ebnmf_tib, idx, glasbey[c(2:11)], "EBNMF (40 components)") +
  theme(axis.text.x = element_blank()) +
  guides(
    fill = guide_legend(ncol = 2)
  )
nmf_p <- do_plot(nmf_tib, idx, glasbey[34:50], "NMF (30 components)") +
  theme(axis.text.x = element_blank()) +
  guides(
    fill = guide_legend(ncol = 2)
  )
nmf40_p <- do_plot(nmf40_tib, idx, glasbey[34:60], "NMF (40 components)") +
  theme(axis.text.x = element_blank()) +
  guides(
    fill = guide_legend(ncol = 2)
  )
tm_p <- do_plot(ebnmf_tib, idx, glasbey[61:80], "Topic Model (30 components)") +
  theme(axis.text.x = element_blank()) +
  guides(
    fill = guide_legend(ncol = 2)
  )
tm40_p <- do_plot(ebnmf_tib, idx, glasbey[61:80], "Topic Model (40 components)") +
  guides(
    fill = guide_legend(ncol = 2)
  )
bas40_p <- plot_grid(ebnmf_p, ebnmf_p, nmf_p, nmf40_p, tm_p, tm40_p,
                     nrow = 6, ncol = 1,
                     rel_heights = c(3, 3, 3, 3, 3, 3.9))
save_plot("./output/plots/montoro_bas40.pdf", bas40_p, base_height = 12.25, base_width = 8)


# Gene plots. -----

LL <- scale_LL(ebnmf$fit$F_pm, ebnmf$fit$L_pm)
exprmean <- log10(readRDS("data/montoro-mean-expr.rds")$var.gene.mean.expr)

k <- 19
ion_tib <- tibble(
  pm = LL[, k],
  SYMBOL = rownames(LL),
  exprmean = exprmean,
  max_other = apply(LL[, -k], 1, max),
  uniqueness = pmax(pmin(pm / max_other, 10), 1)
) |>
  mutate(SYMBOL = ifelse(pm > 0.6 & pm - exprmean > 2.7, SYMBOL, ""))
ggplot(ion_tib, aes(x = exprmean, y = pm, label = SYMBOL, color = uniqueness)) +
  geom_point() +
  geom_text_repel(color = "darkgray",size = 2.25, fontface = "italic",
                  segment.color = "darkgray", segment.size = 0.25,
                  min.segment.length = 0, na.rm = TRUE, max.overlaps = 20) +
  scale_color_gradient(low = "dodgerblue", high = "darkred", trans = "log10") +
  theme_minimal() +
  labs(
    x = "Mean Expression (log10)",
    y = "Gene Score (posterior mean)",
    color = "Uniqueness"
  ) +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    legend.text = element_text(size = 8),
    axis.title = element_text(size = 10, margin = margin(2, 2, 2, 2))
  )
ggsave("output/plots/montoro_ionocyte.pdf", width = 8, height = 4)


k1 <- 15
k2 <- 17
tuft_tib <- tibble(
  cell_type = c(rep("Tuft-1", nrow(LL)), rep("Tuft-2", nrow(LL))),
  pm = c(LL[, k1], LL[, k2]),
  SYMBOL = rep(rownames(LL), times = 2),
  exprmean = rep(exprmean, times = 2),
  max_other = c(apply(LL[, -k1], 1, max), apply(LL[, -k2], 1, max)),
  uniqueness = pmax(pmin(pm / max_other, 10), 1)
) |>
  mutate(SYMBOL = ifelse(pm > 0.85 &
                           (uniqueness > 4 | pm - exprmean > 3),
                         SYMBOL, ""))
ggplot(tuft_tib, aes(x = exprmean, y = pm, label = SYMBOL, color = uniqueness)) +
  geom_point() +
  geom_text_repel(color = "darkgray",size = 2.25, fontface = "italic",
                  segment.color = "darkgray", segment.size = 0.25,
                  min.segment.length = 0, na.rm = TRUE, max.overlaps = 20) +
  scale_color_gradient(low = "dodgerblue", high = "darkred", trans = "log10") +
  theme_minimal() +
  labs(
    x = "Mean Expression (log10)",
    y = "Gene Score (posterior mean)",
    color = "Uniqueness"
  ) +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    legend.text = element_text(size = 8),
    axis.title = element_text(size = 10, margin = margin(2, 2, 2, 2))
  ) +
  facet_wrap(~cell_type, nrow = 1)
ggsave("output/plots/montoro_tuft.pdf", width = 8, height = 4)


k1 <- 24
k2 <- 13
k3 <- 28
goblet_tib <- tibble(
  cell_type = c(rep("Goblet-1 shared", nrow(LL)),
                rep("Goblet-1 larger subset", nrow(LL)),
                rep("Goblet-1 smaller subset", nrow(LL))),
  pm = c(LL[, k1], LL[, k2], LL[, k3]),
  SYMBOL = rep(rownames(LL), times = 3),
  exprmean = rep(exprmean, times = 3),
  max_other = c(apply(LL[, -k1], 1, max),
                apply(LL[, -k2], 1, max),
                apply(LL[, -k3], 1, max)),
  uniqueness = pmax(pmin(pm / max_other, 10), 1)
) |>
  mutate(SYMBOL = ifelse(pm > 1.8 |
                           (uniqueness > 2 & pm > 0.25),
                         SYMBOL, "")) |>
  mutate(cell_type = factor(cell_type, levels = unique(cell_type)))
ggplot(goblet_tib, aes(x = exprmean, y = pm, label = SYMBOL, color = uniqueness)) +
  geom_point() +
  geom_text_repel(color = "darkgray",size = 2.25, fontface = "italic",
                  segment.color = "darkgray", segment.size = 0.25,
                  min.segment.length = 0, na.rm = TRUE, max.overlaps = 20) +
  scale_color_gradient(low = "dodgerblue", high = "darkred", trans = "log10") +
  theme_minimal() +
  labs(
    x = "Mean Expression (log10)",
    y = "Gene Score (posterior mean)",
    color = "Uniqueness"
  ) +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    legend.text = element_text(size = 8),
    axis.title = element_text(size = 10, margin = margin(2, 2, 2, 2))
  ) +
  facet_wrap(~cell_type, nrow = 1)
ggsave("output/plots/montoro_goblet1.pdf", width = 8, height = 4)

k <- 29
gob2_tib <- tibble(
  pm = LL[, k],
  SYMBOL = rownames(LL),
  exprmean = exprmean,
  max_other = apply(LL[, -k], 1, max),
  uniqueness = pmax(pmin(pm / max_other, 10), 1)
) |>
  mutate(SYMBOL = ifelse(pm > 1.8 |
                           (uniqueness > 2 & pm > 0.25),
                         SYMBOL, ""))
ggplot(gob2_tib, aes(x = exprmean, y = pm, label = SYMBOL, color = uniqueness)) +
  geom_point() +
  geom_text_repel(color = "darkgray",size = 2.25, fontface = "italic",
                  segment.color = "darkgray", segment.size = 0.25,
                  min.segment.length = 0, na.rm = TRUE, max.overlaps = 20) +
  scale_color_gradient(low = "dodgerblue", high = "darkred", trans = "log10") +
  theme_minimal() +
  labs(
    x = "Mean Expression (log10)",
    y = "Gene Score (posterior mean)",
    color = "Uniqueness"
  ) +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    legend.text = element_text(size = 8),
    axis.title = element_text(size = 10, margin = margin(2, 2, 2, 2))
  )
ggsave("output/plots/montoro_goblet2.pdf", width = 8, height = 4)


# Top genes (static tables). -----

GSEA <- readRDS("output/montoro-ebnmf-gsea.rds")
keep_res <- GSEA |>
  filter(enrichmentScore > 0) |>
  group_by(k) |>
  filter(FDR < 0.05 | rank(FDR) == min(rank(FDR))) |>
  ungroup()

genesets <- keep_res |>
  mutate(disp_text = paste0(
    description, " (", ifelse(FDR < 0.01, "<0.01", format(round(FDR, 2), nsmall = 2)), ")"
  )) |>
  group_by(k) |>
  summarize(geneSets = paste(disp_text, collapse = ", ")) |>
  mutate(Component = paste0("k", k), .before = k) |>
  select(-k)

LL <- scale_LL(ebnmf$fit$F_pm, ebnmf$fit$L_pm)
colnames(LL) <- paste0("k", 1:ncol(LL))
LL_tib <- as_tibble(LL) |>
  mutate(Gene = gene_names) |>
  pivot_longer(-Gene, names_to = "Component", values_to = "Loading") |>
  mutate(Component = factor(Component, levels = colnames(LL))) |>
  group_by(Gene) |>
  mutate(Uniqueness = Loading / sort(Loading, decreasing = TRUE)[2])

top_unique_genes <- LL_tib |>
  filter(Uniqueness > 2) |>
  group_by(Component) |>
  slice_max(Loading, n = 25) |>
  summarize(TopGenes = paste(Gene, collapse = ", "))

res_tbl <- genesets |>
  left_join(top_unique_genes, by = "Component")

gt_res <- gt(res_tbl) |>
  cols_label(
    geneSets = "Gene Sets (FDR)",
    TopGenes = "Top Unique Genes"
  ) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) |>
  tab_options(table.font.size = px(10))

gtsave(gt_res, "./output/montoro_ebnmf_gene_sets.pdf")
