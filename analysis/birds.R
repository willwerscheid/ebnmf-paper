library(ggplot2)
library(cowplot)
library(tibble)
library(dplyr)
library(tidyr)
source("code/birds_functions.R")

out <- read.csv("data/birds.csv",comment.char = "#",header = FALSE)
out <- as.matrix(out)
d   <- 15
n   <- nrow(out)/d
limbs <- matrix(0,d*d,n)
rows <- 1:d
for (i in 1:n) {
  limbs[,i] <- out[rows,]
  rows <- rows + d
}

# This takes the m x k matrix of feature weights and plots them as n x
# n images, where m = n^2. Note that nrow*ncol >= k is expected.
plot_images <- function (W, n = 15, font_size = 9, nrow = 3, ncol = 3) {
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

# Plot the features.
print(plot_images(limbs))

# Generate all the birds and plot them.
d     <- sqrt(nrow(limbs))
dat   <- expand.grid(left_wing  = c(0,3,5),
                     right_wing = c(0,2,4),
                     left_foot  = c(0,6,8),
                     right_foot = c(0,7,9))
n     <- nrow(dat)
limbs <- cbind(0,limbs)
dat   <- dat + 1
birds <- matrix(0,d*d,n)
torso <- 2
for (i in 1:n) {
  l1 <- dat[i,"left_wing"]
  l2 <- dat[i,"right_wing"]
  l3 <- dat[i,"left_foot"]
  l4 <- dat[i,"right_foot"]
  birds[,i] <- limbs[,torso] + limbs[,l1] + limbs[,l2] +
               limbs[,l3] + limbs[,l4]
}
    
print(plot_images(birds,nrow = 7,ncol = 12))
