# Let 0 < u < 1 be a noise parameter. We simulate values for blank
# pixels (zeroes) from Beta(u,1), and solid pixels (nonzeros) from
# Beta(1,u). The larger the value of u, the noisier the images.
generate_noisy_swimmer_data <- function (dat, u) {
  s1  <- ifelse(dat == 0,u,1)
  s2  <- ifelse(dat == 0,1,u)
  out <- rbeta(length(dat),s1,s2)
  out <- matrix(out,nrow = nrow(dat))
  return(out)
}

# TO DO: Explain what this function does, and how to use it.
frobenius_norm <- function (X, W, H)
  sum((X - W %*% H)^2)

# TO DO: Explain what this function does, and how to use it.
frobenius_norm_nmf <- function (X, fit)
  frobenius_norm(X,fit$W,fit$H)

# TO DO: Explain what this function does, and how to use it.
frobenius_norm_flash <- function (X, fit)
  frobenius_norm(X,fit$L_pm,t(fit$F_pm))

# This takes the m x k matrix of feature weights and plots them as n x
# n images, where m = n^2. Note that nrow*ncol >= k is expected.
plot_images <- function (W, n = 32, font_size = 9, nrow = 4, ncol = 4) {
  k <- ncol(W)
  colnames(W) <- paste0("V",1:k)
  dat <- as_tibble(W)
  dat <- mutate(dat,
                row = rep(1:n,each = n),
                col = rep(1:n,times = n)) 
  dat <- pivot_longer(dat,cols = -c(row,col),names_to = "k",
                      values_to = "loading",names_prefix = "V",
                      names_transform = as.numeric)
  dat$k <- factor(dat$k,1:k)
  return(ggplot(dat,aes(x = row,y = col,fill = loading)) +
    geom_tile() +
    scale_y_reverse() +
    scale_fill_gradient(low = "white",high = "black") +
    facet_wrap(~k,nrow = nrow,ncol = ncol) +
    guides(fill = "none",x = "none",y = "none") +
    labs(x = "",y = "") +
    theme_cowplot(font_size = font_size) +
    theme(strip.background = element_blank(),
          panel.border = element_rect(color = "black",size = 0.5)))
}
