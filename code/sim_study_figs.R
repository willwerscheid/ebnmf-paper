source("./code/sim_functions.R")

ss1_K3 <- readRDS("./output/ss1_K3.rds")
ss1_K4 <- readRDS("./output/ss1_K4.rds")
ss1_K6 <- readRDS("./output/ss1_K6.rds")

ss2_K3 <- readRDS("./output/ss2_K3.rds")
ss2_K4 <- readRDS("./output/ss2_K4.rds")
ss2_K6 <- readRDS("./output/ss2_K6.rds")
ss2_K8 <- readRDS("./output/ss2_K8.rds")


### Plot functions -----

filter_metric <- function(df, metric) {
  df <- df |>
    mutate(method = factor(method, levels = unique(df$method)),
           submethod = factor(submethod, levels = unique(df$submethod))) |>
    group_by(K, method, submethod, varied_n, seed) |>
    summarize(metric_val = max(metric_val * (metric_type %in% metric))) |>
    mutate(metric_val = coalesce(metric_val, 0)) |>
    mutate(varied_n = factor(varied_n)) |>
    group_by(K, method, submethod, varied_n) |>
    summarize(metric_val = mean(metric_val)) |>
    ungroup()
  return(df)
}

make_plots <- function(df, limits, transform, breaks,
                       colors = c("red4", "dodgerblue"),
                       margins = c(38, 1, 32, 40, 40),
                       fill_legend = "Cosine sim.") {
  all_plots <- list()
  methods <- c("NMF", "EBNMF", "KimPark", "RcppML", "Hoyer")

  for (i in 1:length(methods)) {
    next_df <- plot_df |>
      filter(method == methods[[i]]) |>
      mutate(submethod = fct_rev(fct_drop(submethod))) |>
      mutate(metric_val = pmax(pmin(metric_val, limits[2]), limits[1]))
    all_plots[[i]] <- ggplot(next_df, aes(x = varied_n, y = submethod, fill = metric_val)) +
      geom_tile() +
      scale_fill_gradient(low = colors[1], high = colors[2],
                          limits = limits, trans = transform, breaks = breaks) +
      theme_minimal() +
      labs(x = "", y = methods[[i]], fill = fill_legend) +
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
    } else {
      all_plots[[i]] <- all_plots[[i]] +
        theme(axis.text.x = element_text(angle = 45))
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
  filter(!str_detect(submethod, "colwise")) |>
  filter(!(method == "RcppML" & submethod %in% c(0.6, 0.7))) |>
  mutate(submethod = str_remove(submethod, ",.*var"))

plot_df <- filter_metric(df, "LLcosine3")
all_plots <- make_plots(plot_df, c(.82, .9995), transform_logit(), c(.9, .99, .999))

p1 <- plot_grid(plotlist = all_plots$plots, ncol = 1, rel_heights = 3 + c(1.6, 0.8, 8.8, 7, 11))
p2 <- plot_grid(p1, all_plots$legend, nrow = 1, rel_widths = c(9, 1))
p3 <- add_sub(p2, label = "Size of rare population",
              hjust = 0.4, vjust = -0.8, size = 10)
plot(p3)
save_plot("./output/plots/ss1_LL3_allK.pdf", p3, base_width = 9, base_height = 6.5)

## 2. Admixture detection.

df <- ss2_K3 |> mutate(K = "K = 3") |>
  bind_rows(ss2_K4 |> mutate(K = "K = 4")) |>
  bind_rows(ss2_K6 |> mutate(K = "K = 6")) |>
  bind_rows(ss2_K8 |> mutate(K = "K = 8")) |>
  filter(!str_detect(submethod, "colwise")) |>
  filter(!(method == "RcppML" & submethod %in% c(0.6, 0.7))) |>
  mutate(submethod = str_remove(submethod, ",.*var")) |>
  mutate(varied_n = varied_n / 500)

plot_df <- filter_metric(df, c("LLcosine1")) |>
  group_by(K, method, submethod, varied_n) |>
  summarize(metric_val = mean(metric_val)) |>
  ungroup()
all_plots <- make_plots(plot_df, c(.82, .9995), transform_logit(), c(.9, .99, .999))

p1 <- plot_grid(plotlist = all_plots$plots, ncol = 1, rel_heights = 3 + c(1.6, 0.8, 8.8, 7, 11))
p2 <- plot_grid(p1, all_plots$legend, nrow = 1, rel_widths = c(9, 1))
p3 <- add_sub(p2, label = "Proportion of 'pure' individuals",
              hjust = 0.4, vjust = -0.8, size = 10)
plot(p3)
save_plot("./output/plots/ss2_LLall_allK.pdf", p3, base_width = 9, base_height = 6.5)

## 3. Unnecessary factors.

df <- ss2_K6 |> mutate(K = "K = 6") |>
  bind_rows(ss2_K8 |> mutate(K = "K = 8")) |>
  filter(!str_detect(submethod, "alt")) |>
  filter(method %in% c("NMF", "EBNMF") |
           (method == "KimPark" & submethod %in% c("30", "100", "300")) |
           (method == "RcppML" & submethod %in% c("0.01", "0.05", "0.1")) |
           (method == "Hoyer" & submethod %in% c("0.01", "0.05"))) |>
  mutate(varied_n = varied_n / 500)

plot_df <- df |>
  mutate(method = factor(method, levels = unique(df$method)),
         submethod = factor(submethod, levels = unique(df$submethod))) |>
  filter(!(method == "RcppML" & submethod == "0.01")) |>
  group_by(K, method, submethod, varied_n, seed) |>
  filter(str_detect(metric_type, "Scale|elapsed")) |>
  mutate(metric_val = ifelse(str_detect(metric_type, "elapsed"), 0, metric_val)) |>
  summarize(metric_val = (1 - sum(metric_val)) / sum(metric_val)) |>
  mutate(varied_n = factor(varied_n)) |>
  group_by(K, method, submethod, varied_n) |>
  summarize(metric_val = mean(metric_val)) |>
  ungroup()

all_plots <- make_plots(plot_df, c(1, 100), transform_log10(), c(1, 10, 100),
                        margins = c(87, 1, 92, 90, 90),
                        fill_legend = "SNR")

p1 <- plot_grid(plotlist = all_plots$plots, ncol = 1, rel_heights = 3 + c(1.6, 3.5, 2.5, 1, 4))
p2 <- plot_grid(p1, all_plots$legend, nrow = 1, rel_widths = c(9, 1))
p3 <- add_sub(p2, label = "Proportion of 'pure' individuals",
              hjust = 0.4, vjust = -0.8, size = 10)
plot(p3)
save_plot("./output/plots/ss2_SNR.pdf", p3, base_width = 9, base_height = 4)
