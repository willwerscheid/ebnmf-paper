library(tools)
library(flashier)
res <- readRDS("output/montoro-ebmf-12.rds")
fl_k12 <- list(timing = res$t,
               fl_ldf = ldf(res$fit,type = "i"))
fl_k40 <- NULL
save(list = c("fl_k12","fl_k40"),
     file = "montoro_ebmf_fits.rds")
resaveRdaFiles("montoro_ebmf_fits.rds")
