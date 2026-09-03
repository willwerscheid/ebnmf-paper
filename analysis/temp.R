celltypes <- meta_data$celltype
celltypes <- as.character(meta_data$celltype)
celltypes[celltypes == "Club (hillock-associated)"] <- "Club (hillock)"
celltypes[celltypes == "Goblet.1"] <- "Goblet"
celltypes[celltypes == "Goblet.2"] <- "Goblet"
celltypes[celltypes == "Goblet.progenitor"] <- "Goblet"
celltypes[celltypes == "Tuft.1"] <- "Tuft"
celltypes[celltypes == "Tuft.2"] <- "Tuft"
celltypes[celltypes == "Tuft.progenitor"] <- "Tuft"
F <- fl_k30$ldf$F
F <- scale_cols(F)
colnames(F) <- paste0("k",1:30)
rows <- which(meta_data$timepoint == "Tp60")
F <- F[rows,]
celltypes <- celltypes[rows]
round(sort(drop(cor(celltypes == "Ionocyte",F))),digits = 3)
