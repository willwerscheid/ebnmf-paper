source("./code/sim_functions.R")

ss1_K3 <- readRDS("./output/ss1_K3.rds")
ss1_K4 <- readRDS("./output/ss1_K4.rds")
ss1_K6 <- readRDS("./output/ss1_K6.rds")
ss1_K8 <- readRDS("./output/ss1_K8.rds")

ss1_K3_cw <- readRDS("./output/ss1_K3_colwise.rds")
ss1_K4_cw <- readRDS("./output/ss1_K4_colwise.rds")
ss1_K6_cw <- readRDS("./output/ss1_K6_colwise.rds")
ss1_K8_cw <- readRDS("./output/ss1_K8_colwise.rds")

ss2_K3 <- readRDS("./output/ss2_K3.rds")
ss2_K4 <- readRDS("./output/ss2_K4.rds")
ss2_K6 <- readRDS("./output/ss2_K6.rds")
ss2_K8 <- readRDS("./output/ss2_K8.rds")

ss2_K3_cw <- readRDS("./output/ss2_K3_colwise.rds")
ss2_K4_cw <- readRDS("./output/ss2_K4_colwise.rds")
ss2_K6_cw <- readRDS("./output/ss2_K6_colwise.rds")
ss2_K8_cw <- readRDS("./output/ss2_K8_colwise.rds")


### Plot functions -----

filter_metric <- function(df, metric) {
  df <- df |>
    mutate(method = factor(method, levels = unique(df$method)),
           submethod = factor(submethod, levels = unique(df$submethod))) |>
    group_by(K, method, submethod, varied_n, seed) |>
    summarize(metric_val = max(metric_val * (metric_type %in% metric), na.rm = TRUE)) |>
    mutate(metric_val = coalesce(metric_val, 0)) |>
    mutate(varied_n = factor(varied_n)) |>
    group_by(K, method, submethod, varied_n) |>
    summarize(metric_val = mean(metric_val)) |>
    ungroup()
  return(df)
}

filter_scale <- function(df) {
  df |>
    mutate(method = factor(method, levels = unique(df$method)),
           submethod = factor(submethod, levels = unique(df$submethod))) |>
    group_by(K, method, submethod, varied_n, seed) |>
    filter(str_detect(metric_type, "Scale|elapsed")) |>
    filter(metric_type != "Scale4") |>
    mutate(metric_val = ifelse(str_detect(metric_type, "elapsed"), 0, metric_val)) |>
    summarize(metric_val = (1 - sum(metric_val)) / sum(metric_val)) |>
    mutate(varied_n = factor(varied_n)) |>
    group_by(K, method, submethod, varied_n) |>
    summarize(metric_val = mean(metric_val)) |>
    ungroup()
}

make_plots <- function(df, limits, transform, breaks, xlab,
                       colors = c("red4", "dodgerblue", "white"),
                       margins = c(38, 1, 32, 40, 40),
                       fill_legend = "Cosine sim.",
                       ylabs = TRUE,
                       gradient2 = FALSE) {
  all_plots <- list()
  methods <- c("NMF", "EBNMF", "KimPark", "RcppML", "Hoyer")

  for (i in 1:length(methods)) {
    next_df <- df |>
      filter(method == methods[[i]]) |>
      mutate(submethod = fct_rev(fct_drop(submethod))) |>
      mutate(metric_val = pmax(pmin(metric_val, limits[2]), limits[1]))
    if (gradient2) {
      all_plots[[i]] <- ggplot(next_df, aes(x = varied_n, y = submethod, fill = metric_val)) +
        geom_tile() +
        scale_fill_gradient2(low = colors[1], high = colors[2], mid = colors[3], midpoint = 0,
                             limits = limits, trans = transform, breaks = breaks) +
        theme_minimal() +
        labs(x = "", y = methods[[i]], fill = fill_legend) +
        theme(plot.margin = margin(l = margins[i]),
              legend.position = "bottom") +
        facet_wrap(~K, nrow = 1)
    } else {
      all_plots[[i]] <- ggplot(next_df, aes(x = varied_n, y = submethod, fill = metric_val)) +
        geom_tile() +
        scale_fill_gradient(low = colors[1], high = colors[2],
                            limits = limits, trans = transform, breaks = breaks) +
        theme_minimal() +
        labs(x = "", y = methods[[i]], fill = fill_legend) +
        theme(plot.margin = margin(l = margins[i]),
              legend.position = "bottom") +
        facet_wrap(~K, nrow = 1)
    }
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
    if (!ylabs) {
      all_plots[[i]] <- all_plots[[i]] +
        labs(y = "") +
        theme(axis.text.y = element_blank())
    }
    if (i < length(methods)) {
      all_plots[[i]] <- all_plots[[i]] +
        theme(axis.text.x = element_blank())
    } else {
      all_plots[[i]] <- all_plots[[i]] +
        theme(axis.text.x = element_text(angle = 45)) +
        labs(x = xlab)
    }
  }

  legend <- get_legend(all_plots[[1]])
  for (i in 1:length(methods)) {
    all_plots[[i]] <- all_plots[[i]] +
      theme(legend.position = "none")
  }

  return(list(plots = all_plots, legend = legend))
}


### Plots -----


## 1. Rare population detection.

df <- ss1_K3 |> mutate(K = "K = 3") |>
  bind_rows(ss1_K4 |> mutate(K = "K = 4")) |>
  bind_rows(ss1_K6 |> mutate(K = "K = 6")) |>
  bind_rows(ss1_K8 |> mutate(K = "K = 8")) |>
  filter(!(method == "RcppML" & submethod %in% c(0.6, 0.7)))
plot_df <- filter_metric(df, "LLcosine3")
all_plots <- make_plots(plot_df, c(.1, .999), transform_logit(), c(.5, .9, .99),
                        xlab = "Size of small cluster",
                        colors = c("red4", "lightyellow"),
                        margins = c(87, 2, 83, 90, 90))

p1 <- plot_grid(plotlist = all_plots$plots, ncol = 1, rel_heights = 3 + c(1, 3.2, 8.8, 7, 11)) #c(1.6, 0.8, 8.8, 7, 11))
p2 <- plot_grid(p1, all_plots$legend, ncol = 1, rel_heights = c(9, 1))
plot(p2)

# df2 <- df |> filter(K %in% c("K = 6", "K = 8"))
# plot_df2 <- filter_scale(df2)
# all_plots2 <- make_plots(plot_df2, c(1, 100), transform_log10(), c(1, 10, 100),
#                          xlab = "Size of small cluster",
#                          margins = rep(30, 5),
#                          fill_legend = "SNR",
#                          ylabs = FALSE)

df2 <- df
plot_df2 <- filter_metric(df2, "SparsenessRatio")
all_plots2 <- make_plots(plot_df2, c(.33, 1.5), transform_log10(), c(.5, 1),
                         xlab = "Size of small cluster",
                         margins = rep(30, 5),
                         fill_legend = "Sparseness (ratio to true)",
                         colors = c("red4", "blue4", "lightyellow"),
                         ylabs = FALSE, gradient2 = TRUE)

p3 <- plot_grid(plotlist = all_plots2$plots, ncol = 1, rel_heights = 3 + c(1, 3.2, 8.8, 7, 11)) #rel_heights = 3 + c(1.6, 3.5, 2.5, 1, 4))
p4 <- plot_grid(p3, all_plots2$legend, ncol = 1, rel_heights = c(9, 1))
plot(p4)

p <- plot_grid(p2, p4, nrow = 1, rel_widths = c(4, 3.1),
               labels = letters[1:2], label_x = 0.1)

save_plot("./output/plots/ss1_LL3_allK.pdf", p, base_width = 11, base_height = 8)

## 1b. Rare population detection (column-wise variance structure).

df <- ss1_K3_cw |> mutate(K = "K = 3") |>
  bind_rows(ss1_K4_cw |> mutate(K = "K = 4")) |>
  bind_rows(ss1_K6_cw |> mutate(K = "K = 6")) |>
  bind_rows(ss1_K8_cw |> mutate(K = "K = 8")) |>
  filter(!(method == "RcppML" & submethod %in% c(0.6, 0.7))) |>
  filter(varied_n >= 5)
plot_df <- filter_metric(df, "LLcosine3")
all_plots <- make_plots(plot_df, c(.1, .999), transform_logit(), c(.5, .9, .99),
                        xlab = "Size of small cluster",
                        colors = c("red4", "lightyellow"),
                        margins = c(87, 2, 83, 90, 90))

p1 <- plot_grid(plotlist = all_plots$plots, ncol = 1, rel_heights = 3 + c(1, 3.2, 8.8, 7, 11)) #c(1.6, 0.8, 8.8, 7, 11))
p2 <- plot_grid(p1, all_plots$legend, ncol = 1, rel_heights = c(9, 1))
plot(p2)

df2 <- df
plot_df2 <- filter_metric(df2, "SparsenessRatio")
all_plots2 <- make_plots(plot_df2, c(.33, 1.5), transform_log10(), c(.5, 1),
                         xlab = "Size of small cluster",
                         margins = rep(30, 5),
                         fill_legend = "Sparseness (ratio to true)",
                         colors = c("red4", "blue4", "lightyellow"),
                         ylabs = FALSE, gradient2 = TRUE)

p3 <- plot_grid(plotlist = all_plots2$plots, ncol = 1, rel_heights = 3 + c(1, 3.2, 8.8, 7, 11)) #rel_heights = 3 + c(1.6, 3.5, 2.5, 1, 4))
p4 <- plot_grid(p3, all_plots2$legend, ncol = 1, rel_heights = c(9, 1))
plot(p4)

p <- plot_grid(p2, p4, nrow = 1, rel_widths = c(4, 3.1),
               labels = letters[1:2], label_x = 0.1)

save_plot("./output/plots/ss1_cw_LL3_allK.pdf", p, base_width = 11, base_height = 8)


## 2. Admixture detection.

df <- ss2_K3 |> mutate(K = "K = 3") |>
  bind_rows(ss2_K4 |> mutate(K = "K = 4")) |>
  bind_rows(ss2_K6 |> mutate(K = "K = 6")) |>
  bind_rows(ss2_K8 |> mutate(K = "K = 8")) |>
  filter(!(method == "RcppML" & submethod %in% c(0.6, 0.7))) |>
  mutate(varied_n = varied_n / 500)

plot_df <- filter_metric(df, c("LLcosine1")) |>
  group_by(K, method, submethod, varied_n) |>
  summarize(metric_val = mean(metric_val)) |>
  ungroup()
all_plots <- make_plots(plot_df, c(.1, .999), transform_logit(), c(.5, .9, .99),
                        xlab = "Proportion of 'pure' observations",
                        colors = c("red4", "lightyellow"),
                        margins = c(87, 2, 83, 90, 90))

p1 <- plot_grid(plotlist = all_plots$plots, ncol = 1, rel_heights = 3 + c(1, 3.2, 8.8, 7, 11)) #c(1.6, 0.8, 8.8, 7, 11))
p2 <- plot_grid(p1, all_plots$legend, ncol = 1, rel_heights = c(9, 1))
plot(p2)

df2 <- df
plot_df2 <- filter_metric(df2, "SparsenessRatio")
all_plots2 <- make_plots(plot_df2, c(.33, 1.5), transform_log10(), c(.5, 1),
                         xlab = "Proportion of 'pure' observations",
                         margins = rep(30, 5),
                         fill_legend = "Sparseness (ratio to true)",
                         colors = c("red4", "blue4", "lightyellow"),
                         ylabs = FALSE, gradient2 = TRUE)

p3 <- plot_grid(plotlist = all_plots2$plots, ncol = 1, rel_heights = 3 + c(1, 3.2, 8.8, 7, 11)) #rel_heights = 3 + c(1.6, 3.5, 2.5, 1, 4))
p4 <- plot_grid(p3, all_plots2$legend, ncol = 1, rel_heights = c(9, 1))
plot(p4)

p <- plot_grid(p2, p4, nrow = 1, rel_widths = c(4, 3.1),
               labels = letters[1:2], label_x = 0.1)

save_plot("./output/plots/ss2_LLall_allK.pdf", p, base_width = 11, base_height = 8)

## 2b. Admixture detection (column-wise variance structure).

df <- ss2_K3_cw |> mutate(K = "K = 3") |>
  bind_rows(ss2_K4_cw |> mutate(K = "K = 4")) |>
  bind_rows(ss2_K6_cw |> mutate(K = "K = 6")) |>
  bind_rows(ss2_K8_cw |> mutate(K = "K = 8")) |>
  filter(!(method == "RcppML" & submethod %in% c(0.6, 0.7))) |>
  mutate(varied_n = varied_n / 500)

plot_df <- filter_metric(df, c("LLcosine1")) |>
  group_by(K, method, submethod, varied_n) |>
  summarize(metric_val = mean(metric_val)) |>
  ungroup()
all_plots <- make_plots(plot_df, c(.1, .999), transform_logit(), c(.5, .9, .99),
                        xlab = "Proportion of 'pure' observations",
                        colors = c("red4", "lightyellow"),
                        margins = c(87, 2, 83, 90, 90))

p1 <- plot_grid(plotlist = all_plots$plots, ncol = 1, rel_heights = 3 + c(1, 3.2, 8.8, 7, 11)) #c(1.6, 0.8, 8.8, 7, 11))
p2 <- plot_grid(p1, all_plots$legend, ncol = 1, rel_heights = c(9, 1))
plot(p2)

df2 <- df
plot_df2 <- filter_metric(df2, "SparsenessRatio")
all_plots2 <- make_plots(plot_df2, c(.33, 1.5), transform_log10(), c(.5, 1),
                         xlab = "Proportion of 'pure' observations",
                         margins = rep(30, 5),
                         fill_legend = "Sparseness (ratio to true)",
                         colors = c("red4", "blue4", "lightyellow"),
                         ylabs = FALSE, gradient2 = TRUE)

p3 <- plot_grid(plotlist = all_plots2$plots, ncol = 1, rel_heights = 3 + c(1, 3.2, 8.8, 7, 11)) #rel_heights = 3 + c(1.6, 3.5, 2.5, 1, 4))
p4 <- plot_grid(p3, all_plots2$legend, ncol = 1, rel_heights = c(9, 1))
plot(p4)

p <- plot_grid(p2, p4, nrow = 1, rel_widths = c(4, 3.1),
               labels = letters[1:2], label_x = 0.1)

save_plot("./output/plots/ss2_cw_LLall_allK.pdf", p, base_width = 11, base_height = 8)


## 4. Timing comparisons.

t_df <- ss2_K3 |> mutate(K = "K = 3") |>
  bind_rows(ss2_K4 |> mutate(K = "K = 4")) |>
  bind_rows(ss2_K6 |> mutate(K = "K = 6")) |>
  bind_rows(ss2_K8 |> mutate(K = "K = 8")) |>
  filter(!(method == "RcppML" & submethod %in% c(0.6, 0.7)))

plot_df <- filter_metric(t_df, c("t_elapsed")) |>
  mutate(method = as.character(method), submethod = as.character(submethod)) |>
  mutate(method = ifelse(method == "NMF", "RcppML", method)) |>
  mutate(method = ifelse(method == "RcppML", "RcppML (includes vanilla NMF)", method)) |>
  mutate(method = ifelse(method == "EBNMF", paste0(method, " (", submethod, ")"), method)) |>
  mutate(method = factor(method, levels = c("EBNMF (NMF init, const var)",
                                            "EBNMF (NMF init, colwise var)",
                                            "EBNMF (greedy/backfit, const var)",
                                            "EBNMF (greedy/backfit, colwise var)",
                                            "RcppML (includes vanilla NMF)",
                                            "KimPark",
                                            "Hoyer"))) |>
  ungroup()

ggplot(plot_df, aes(x = K, y = metric_val, color = method)) +
  geom_boxplot(outliers = FALSE) +
  scale_y_log10() +
  scale_color_brewer(palette = "Set2") +
  labs(x = "", y = "Elapsed time (s)", color = "") +
  guides(color = guide_legend(nrow = 4, ncol = 2, position = "bottom", direction = "horizontal")) +
  theme_minimal()

ggsave("./output/plots/ss2_t_elapsed.pdf", width = 7, height = 4.5)
