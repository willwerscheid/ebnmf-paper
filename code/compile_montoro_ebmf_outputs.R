library(tools)
library(flashier)
res <- readRDS("output/montoro-ebmf-12.rds")
fl_k12 <- list(timing = res$t,
               fl_ldf = ldf(res$fit,type = "i"),
               L_ghat = res$fit$L_ghat,
               F_ghat = res$fit$F_ghat,
               residuals_sd = res$fit$residuals_sd)
res <- readRDS("output/montoro-ebmf-40.rds")
fl_k40 <- list(timing = res$t,
               fl_ldf = ldf(res$fit,type = "i"),
               L_ghat = res$fit$L_ghat,
               F_ghat = res$fit$F_ghat,
               residuals_sd = res$fit$residuals_sd)
save(list = c("fl_k12","fl_k40"),
     file = "montoro_ebmf_fits.RData")
resaveRdaFiles("montoro_ebmf_fits.RData")
