# TO DO: Explain here what this script is for, and how to use it.
#
# The CBCL faces data set was downloaded from here:
# http://www.ai.mit.edu/courses/6.899/lectures/faces.tar.gz
#
# In each file, the first line contains the number of examples (6,977
# training, 24,045 test), the second line contains the number of
# dimensions (19x19 = 361). Each additional line consists of a single
# image, histogram equalized and normalized so that all pixel values
# are between 0 and 1, followed by a 1 for a face or a -1 for a
# non-face.
dat <- readLines("svm.train.normgrey")

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
i <- which(image_labels == 1)
image_labels <- image_labels[i]
X <- X[i,]
