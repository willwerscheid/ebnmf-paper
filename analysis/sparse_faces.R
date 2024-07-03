W       <- normalize.cols(nmf$W)
L_flnmf <- ldf(fl_nmf,type = "i")$L
L_fl    <- ldf(fl,type = "i")$L
cor(c(W),c(L_flnmf))
cor(c(W),c(L_fl))

zero    <- 0.001
k_set   <- order(colSums(W > zero))
W       <- W[,k_set]
L_flnmf <- L_flnmf[,k_set]
L_fl    <- L_fl[,k_set]
plot_faces(L_flnmf,title = "flashier with good NMF init")
plot_faces(L_fl,title = "flashier with rough NMF init")

# Compare sparsity of the solutions.
pdat <-
  data.frame(
    nmf       = sort(colMeans(normalize.cols(W) > zero)),
    flashier1 = sort(colMeans(normalize.cols(L_flnmf) > zero)),
    flashier2 = sort(colMeans(normalize.cols(L_fl) > zero)))
p1 <- ggplot(pdat,aes(x = nmf,y = flashier1)) +
  geom_point() +
  geom_abline(intercept = 0,slope = 1,color = "magenta",lty = "dotted") +
  theme_cowplot(font_size = 12)
p2 <- ggplot(pdat,aes(x = nmf,y = flashier2)) +
  geom_point() +
  geom_abline(intercept = 0,slope = 1,color = "magenta",lty = "dotted") +
  theme_cowplot(font_size = 12)
plot_grid(p1,p2)
