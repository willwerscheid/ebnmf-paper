library(tidyverse)

go <- c("geneontology_Biological_Process",
        "geneontology_Cellular_Component",
        "geneontology_Molecular_Function")

run_GSEA <- function(LL, incl_thresh = 0.01) {
  all_res <- tibble()
  for (k in 1:ncol(LL)) {
    dat <- data.frame(
      Gene = names(LL[, k]),
      Loading = LL[, k]
    ) |>
      filter(Loading > max(Loading) * incl_thresh)
    cat("Running GSEA on component", k,
        "with", nrow(dat), "genes...\n")
    wgres <- WebGestaltR::WebGestaltR(
      organism = "mmusculus",
      enrichMethod = "GSEA",
      enrichDatabase = go,
      interestGene = dat, interestGeneType = "genesymbol",
      minNum = 10, maxNum = 500,
      topThr = 10, sigMethod = "top",
      perNum = 100,
      isOutput = FALSE
    )
    wgres <- wgres |>
      mutate(k = k, .before = geneSet)
    all_res <- all_res |>
      bind_rows(wgres)
  }
  return(all_res)
}

ebnmf <- readRDS("output/montoro-ebmf.rds")
ebnmf_LL <- ebnmf$fit$L_pm
ebnmf_GSEA <- run_GSEA(ebnmf_LL)
saveRDS(ebnmf_GSEA, "output/montoro-ebnmf-gsea.rds")
