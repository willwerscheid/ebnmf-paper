# The CBCL faces data set was downloaded from here:
# http://www.ai.mit.edu/courses/6.899/lectures/faces.tar.gz
#
# In each file, the first line contains the number of examples (6,977
# training, 24,045 test), the second line contains the number of
# dimensions (19x19 = 361). Each additional line consists of a single
# image, histogram equalized and normalized so that all pixel values
# are between 0 and 1, followed by a 1 for a face or a -1 for a
# non-face.
library(tools)
source("../code/faces_functions.R")
faces_train <- read_faces_data("svm.train.normgrey")
faces_test <- read_faces_data("svm.test.normgrey")
save(list = c("faces_train","faces_test"),
     file = "faces.RData")
resaveRdaFiles("faces.RData")
