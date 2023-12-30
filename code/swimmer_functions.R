# TO DO: Explain here what this function is for, and how to use it.
plot_images <- function (W) {
  W <- max(W) - W
  colnames(W) <- paste0("V", 1:ncol(W))
  tib <- as_tibble(W)
  tib <- mutate(tib,
    row = rep(1:32,each = 32),
    col = rep(1:32,times = 32)) 
  tib <- pivot_longer(tib,
      cols = -c(row,col),
      names_to = "k",
      values_to = "loading",
      names_prefix = "V",
      names_transform = as.numeric)
  p <- ggplot(tib,aes(x = row,y = col,fill = loading)) +
    geom_tile() +
    scale_y_reverse() +
    scale_fill_gradient(low = "white",high = "black") +
    facet_wrap(~k) +
    guides(fill = "none",x = "none",y = "none") +
    theme_cowplot(font_size = 8)
  return(p)
}
