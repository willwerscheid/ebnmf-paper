library(Matrix)

preprocess <- function(dat, min.nzcts = 10) {
  size.factors <- Matrix::colSums(dat)
  size.factors <- size.factors / mean(size.factors)
  pc <- 1

  # Use liger for most variable gene selection.
  liger.dat <- rliger::createLiger(list(dat = dat))
  liger.dat <- rliger::normalize(liger.dat)
  liger.dat <- rliger::selectGenes(liger.dat)

  gene_cts <- Matrix::rowSums(dat > 0)
  dat <- dat[gene_cts >= min.nzcts, ]

  return(list(
    counts = dat,
    sf = size.factors,
    pc = pc,
    var.genes = intersect(liger.dat@varFeatures, rownames(dat))
  ))
}

zz <- download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE103354&format=file&file=GSE103354%5FPulseSeq%5FUMI%5Fcounts%2Erds%2Egz",
                    destfile = "./data/pulseseq.rds.gz")
R.utils::gunzip("./data/pulseseq.rds.gz")

dat <- readRDS("./data/pulseseq.rds")
file.remove("./data/pulseseq.rds")

source("./code/preprocess.R")
pp.dat <- preprocess(dat)
saveRDS(pp.dat, "./data/montoro.rds")

mean.expr <- rowSums(pp.dat$counts) / ncol(pp.dat$counts)
var.gene.mean.expr <- mean.expr[pp.dat$var.genes]
saveRDS(
  list(mean.expr = mean.expr, var.gene.mean.expr = var.gene.mean.expr),
  "./data/montoro-mean-expr.rds"
)
