# Ad hoc method to extract the true factors from the (non-noisy)
# swimmer data set.
#
# Note that Y should be loaded into R by taking the following steps:
#
#   Y <- readMat("data/swimmer.mat")$Y
#   Y <- apply(Y,3,as.vector)
#   Y <- Y - 1
#
extract_true_swimmer_factors <- function (Y) {
  Y <- Y + 0.0001  
  fit <- flash(Y,ebnm_fn = ebnm_point_exponential,
               backfit = FALSE,var_type = 0)
  L <- ldf(fit)$L
  a <- as.numeric(rowSums(L[,-1]) < 0.01)
  L[,1] <- a * L[,1]
  L <- (L > 0.1)
  mode(L) <- "numeric"
  return(L)
}

# Let 0 < u < 1 be a noise parameter. We simulate values for blank
# pixels (zeroes) from Beta(u,1), and solid pixels (ones) from
# Beta(1,u). The larger the value of u, the noisier the images.
generate_noisy_swimmer_data <- function (dat, u) {
  s1  <- ifelse(dat == 0,u,1)
  s2  <- ifelse(dat == 0,1,u)
  out <- rbeta(length(dat),s1,s2)
  out <- matrix(out,nrow = nrow(dat))
  return(out)
}

# Compute the Frobenius norm objective for data matrix X and matrix
# factorization W*H.
frobenius_norm <- function (X, W, H)
  sum((X - W %*% H)^2)

# Compute the Frobenius norm objective for data matrix X and NMF
# result fit$W * fit$H.
frobenius_norm_nmf <- function (X, fit)
  frobenius_norm(X,fit$W,fit$H)

# Compute the Frobenius norm objective for data matrix X and flashier
# output "fit".
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
