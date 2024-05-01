# This function reads the data from svm.train.normgrey or
# svm.test.normgrey and outputs an n x m matrix, where n is the number
# of faces and m is the number of pixels.
read_faces_data <- function (filename) {

  # Read from the text file.
  dat <- readLines(filename)

  # Get the number of images.
  n <- as.numeric(dat[[1]])

  # Get the number of pixels in each image.
  m <- as.numeric(dat[[2]])

  # Initialize the n x m data matrix and the face/non-face labels.
  X <- matrix(0,n,m)
  image_labels <- rep(0,n)

  # Extract the images.
  dat <- dat[-(1:2)]
  for (i in 1:n) {
    pixels <- as.numeric(unlist(strsplit(trimws(dat[[i]])," ")))
    image_labels[i] <- pixels[m+1]
    X[i,] <- pixels[1:m]
  }

  # Keep only the images that are faces.
  # 1 = face, -1 = non-face.
  return(X[image_labels == 1,])
}
