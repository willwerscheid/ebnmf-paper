# Normalize the columns of W so that the maximum value in each column
# is 1.
normalize.cols <- function (W)
  W <- apply(W,2,function (x) x/max(x))

# Compute the Frobenius norm objective for data matrix X and matrix
# factorization W*H.
frobenius_norm <- function (X, W, H)
  sum((X - W %*% H)^2)

# This function reads the data from svm.train.normgrey or
# svm.test.normgrey and outputs an m x n matrix, where n is the number
# of pixels and m is the number of faces.
read_faces_data <- function (filename) {

  # Read from the text file.
  dat <- readLines(filename)

  # Get the number of images.
  n <- as.numeric(dat[[1]])

  # Get the number of pixels in each image.
  m <- as.numeric(dat[[2]])

  # Initialize the n x m data matrix and the face/non-face labels.
  X <- matrix(0,m,n)
  image_labels <- rep(0,n)

  # Extract the images.
  dat <- dat[-(1:2)]
  for (i in 1:n) {
    pixels <- as.numeric(unlist(strsplit(trimws(dat[[i]])," ")))
    image_labels[i] <- pixels[m+1]
    X[,i] <- pixels[1:m]
  }

  # Keep only the images that are faces.
  # 1 = face, -1 = non-face.
  return(X[,image_labels == 1])
}

# This takes the m x k matrix of feature weights and plots them as n x
# n images, where m = n^2. Note that nrow*ncol >= k is expected.
plot_faces <- function (W, n = 19, font_size = 9, nrow = 5, ncol = 10,
                        title = "") {
  k <- ncol(W)
  W <- normalize.cols(W)
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
    coord_flip() +
    scale_x_reverse() +
    scale_fill_gradient(low = "white",high = "black") +
    facet_wrap(~k,nrow = nrow,ncol = ncol) +
    guides(fill = "none",x = "none",y = "none") +
    labs(x = "",y = "",title = title) +
    theme_cowplot(font_size = font_size) +
    theme(strip.background = element_blank(),
          panel.border = element_rect(color = "black",linewidth = 0.5),
          plot.title = element_text(size = font_size,face = "plain")))
}
