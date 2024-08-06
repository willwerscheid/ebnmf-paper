# Generate a 225 x 10 matrix containing the bird limbs (15 x 15 = 225).
read_bird_limbs <- function (filename) {
  d    <- 15
  dat  <- read.csv(filename,comment.char = "#",header = FALSE)
  dat  <- as.matrix(dat)
  n    <- nrow(dat)/d
  out  <- matrix(0,d*d,n)
  rows <- 1:d
  for (i in 1:n) {
    out[,i] <- dat[rows,]
    rows <- rows + d 
  }
  return(out)
}

# This takes the m x k matrix of feature weights and plots them as n x
# n images, where m = n^2. Note that nrow*ncol >= k is expected.
plot_images <- function (W, n = 15, font_size = 9, nrow = 9, ncol = 9) {
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
          panel.border = element_rect(color = "black",linewidth = 0.5)))
}

# Generate birds from all possible combinations of the limbs; the
# output should be a 225 x 81 matrix.
generate_all_birds <- function (limbs) {
  d     <- sqrt(nrow(limbs))
  dat   <- expand.grid(left_wing  = c(0,3,5),
                       right_wing = c(0,2,4),
                       left_foot  = c(0,6,8),
                       right_foot = c(0,7,9))
  n     <- nrow(dat)
  limbs <- cbind(0,limbs)
  dat   <- dat + 1
  out   <- matrix(0,d*d,n)
  torso <- 2
  for (i in 1:n) {
    l1 <- dat[i,"left_wing"]
    l2 <- dat[i,"right_wing"]
    l3 <- dat[i,"left_foot"]
    l4 <- dat[i,"right_foot"]
    out[,i] <- limbs[,torso] +
               limbs[,l1] + limbs[,l2] + limbs[,l3] + limbs[,l4]
  }
  return(out)
}
