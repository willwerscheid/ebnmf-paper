library(tidyverse)
library(cowplot)
library(ggrepel)
library(ggpubr)
library(reactable)

glasbey <- fastTopics:::glasbey()[1:200]

scale_FF <- function(FF, LL, D = 1) {
  FF_norms <- apply(FF, 2, max)
  FF <- t(t(FF) / FF_norms)
  return(FF)
}
scale_LL <- function(FF, LL, D = 1) {
  FF_norms <- apply(FF, 2, max) * D
  LL <- t(t(LL) * FF_norms)
  return(LL)
}

ebnmf <- readRDS("output/montoro-ebmf.rds")
ebnmf_FF <- scale_FF(ebnmf$fit$F_pm, ebnmf$fit$L_pm)

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
  CellIdx = 1:nrow(ebnmf_FF),
  CellType = sapply(strsplit(rownames(ebnmf_FF), "_"), `[[`, 5),
  TimePoint = sapply(strsplit(rownames(ebnmf_FF), "_"), `[[`, 2),
  Replicate = sapply(strsplit(rownames(ebnmf_FF), "_"), `[[`, 3)
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

get_top_genes <- function(LL, min_loading = 0.5, max_n = 15) {
  colnames(LL) <- paste0("k", 1:ncol(LL))
  LL_tib <- as_tibble(LL) |>
    mutate(Gene = rownames(LL)) |>
    pivot_longer(-Gene, names_to = "Component", values_to = "Loading") |>
    mutate(Component = factor(Component, levels = colnames(LL))) |>
    group_by(Gene) |>
    mutate(Uniqueness = Loading / sort(Loading, decreasing = TRUE)[2])

  top_unique_genes <- LL_tib |>
    filter(Loading > min_loading, Uniqueness > 1) |>
    mutate(Gene = str_remove(Gene, "_.*")) |>
    group_by(Component) |>
    slice_max(Uniqueness, n = max_n, with_ties = FALSE) |>
    summarize(TopGenes = paste(Gene, collapse = ", "))

  top_unique_genes <- top_unique_genes |>
    mutate(ComponentName = paste0(Component, " (", TopGenes, ")"))

  return(top_unique_genes)
}


# Main figures. -----

make_main_fig <- function(LL, FF, cell_idx, p_title, min_loading = 0.5, max_n = 15) {
  top_unique_genes <- get_top_genes(LL, min_loading, max_n)

  ct_tib <- make_tib(FF, cell_idx)
  ct_tib_idx <- do_order(ct_tib)

  plot_df <- ct_tib |>
    left_join(top_unique_genes, by = "Component") |>
    mutate(Component = factor(
      ComponentName, levels = top_unique_genes$ComponentName)
    ) |>
    select(-TopGenes, -ComponentName)

  p <- do_plot(plot_df, ct_tib_idx, glasbey[2:41], p_title, legend_ncol = 1) +
    labs(fill = "Component (Top unique genes)") +
    theme(legend.position = "bottom",
          legend.title.position = "top",
          legend.title = element_text(face = "bold"))

  return(p)
}

set.seed(1)
samp_idx <- external_info |>
  slice_sample(n = 200, by = CellType, replace = TRUE) |>
  pull(CellIdx)

ebnmf_LL <- scale_LL(ebnmf$fit$F_pm, ebnmf$fit$L_pm)
ebnmf_p <- make_main_fig(
  ebnmf_LL, ebnmf_FF, samp_idx, "EBNMF results (K = 30)"
)
save_plot("./output/plots/montoro_ebnmf.pdf",
          ebnmf_p, base_height = 10, base_width = 8)

ebnmf40 <- readRDS("output/montoro-ebmf-40.rds")
ebnmf40_FF <- scale_FF(ebnmf40$fit$F_pm, ebnmf40$fit$L_pm)
ebnmf40_LL <- scale_LL(ebnmf40$fit$F_pm, ebnmf40$fit$L_pm)
ebnmf40_p <- make_main_fig(
  ebnmf40_LL, ebnmf40_FF, samp_idx, "EBNMF results (K = 40)"
)
save_plot("./output/plots/montoro_ebnmf40.pdf",
          ebnmf40_p, base_height = 12, base_width = 8)

nmf <- readRDS("output/montoro-nmf.rds")
nmf_FF <- scale_FF(t(nmf$fit$h), nmf$fit$w, nmf$fit$d)
nmf_LL <- scale_LL(t(nmf$fit$h), nmf$fit$w, nmf$fit$d)
rownames(nmf_LL) <- rownames(ebnmf_LL)
nmf_p <- make_main_fig(
  nmf_LL, nmf_FF, samp_idx, "NMF results (K = 30)", max_n = 14
)
save_plot("./output/plots/montoro_nmf.pdf",
          nmf_p, base_height = 10, base_width = 8)

nmf40 <- readRDS("output/montoro-nmf-40.rds")
nmf40_FF <- scale_FF(t(nmf40$fit$h), nmf40$fit$w, nmf40$fit$d)
nmf40_LL <- scale_LL(t(nmf40$fit$h), nmf40$fit$w, nmf40$fit$d)
rownames(nmf40_LL) <- rownames(ebnmf_LL)
nmf40_p <- make_main_fig(
  nmf40_LL, nmf40_FF, samp_idx, "NMF results (K = 40)", max_n = 14
)
save_plot("./output/plots/montoro_nmf40.pdf",
          nmf40_p, base_height = 12, base_width = 8)

tm <- readRDS("output/montoro-topics.rds")
tm_FF <- scale_FF(tm$fit$F, tm$fit$L)
tm_LL <- scale_LL(tm$fit$F, tm$fit$L)
tm_p <- make_main_fig(
  tm_LL, tm_FF, samp_idx, "Topic model results (K = 30)", max_n = 14
)
save_plot("./output/plots/montoro_tm.pdf",
          tm_p, base_height = 10, base_width = 8)

tm40 <- readRDS("output/montoro-topics-40.rds")
tm40_FF <- scale_FF(tm40$fit$F, tm40$fit$L)
tm40_LL <- scale_LL(tm40$fit$F, tm40$fit$L)
tm40_p <- make_main_fig(
  tm40_LL, tm40_FF, samp_idx, "Topic model results (K = 40)", max_n = 14
)
save_plot("./output/plots/montoro_tm40.pdf",
          tm40_p, base_height = 12, base_width = 8)

combined_p <- plot_grid(
  ebnmf_p + theme(legend.position = "none"),
  ebnmf40_p + theme(legend.position = "none"),
  nmf_p + theme(legend.position = "none"),
  nmf40_p + theme(legend.position = "none"),
  tm_p + theme(legend.position = "none"),
  tm40_p + theme(legend.position = "none"),
  ncol = 1
)
save_plot("./output/plots/montoro_combined.pdf",
          combined_p, base_height = 12, base_width = 8)


# Gene plots. -----

exprmean <- log10(readRDS("data/montoro-mean-expr.rds")$var.gene.mean.expr)

k <- 19
ion_tib <- tibble(
  pm = ebnmf_LL[, k],
  SYMBOL = rownames(ebnmf_LL),
  exprmean = exprmean,
  max_other = apply(ebnmf_LL[, -k], 1, max),
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
  cell_type = c(rep("Tuft-1", nrow(ebnmf_LL)), rep("Tuft-2", nrow(ebnmf_LL))),
  pm = c(ebnmf_LL[, k1], ebnmf_LL[, k2]),
  SYMBOL = rep(rownames(ebnmf_LL), times = 2),
  exprmean = rep(exprmean, times = 2),
  max_other = c(apply(ebnmf_LL[, -k1], 1, max), apply(ebnmf_LL[, -k2], 1, max)),
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
  cell_type = c(rep("Goblet-1 shared", nrow(ebnmf_LL)),
                rep("Goblet-1 larger subset", nrow(ebnmf_LL)),
                rep("Goblet-1 smaller subset", nrow(ebnmf_LL))),
  pm = c(ebnmf_LL[, k1], ebnmf_LL[, k2], ebnmf_LL[, k3]),
  SYMBOL = rep(rownames(ebnmf_LL), times = 3),
  exprmean = rep(exprmean, times = 3),
  max_other = c(apply(ebnmf_LL[, -k1], 1, max),
                apply(ebnmf_LL[, -k2], 1, max),
                apply(ebnmf_LL[, -k3], 1, max)),
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
  pm = ebnmf_LL[, k],
  SYMBOL = rownames(ebnmf_LL),
  exprmean = exprmean,
  max_other = apply(ebnmf_LL[, -k], 1, max),
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

ebnmf_LL <- scale_LL(ebnmf$fit$F_pm, ebnmf$fit$L_pm)
colnames(ebnmf_LL) <- paste0("k", 1:ncol(ebnmf_LL))
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
