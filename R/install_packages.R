if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery", "limma", "org.Hs.eg.db", "WGCNA", "dorothea", "viper",
                       "TCGAbiolinks", "SummarizedExperiment", "PharmacoGx", "Biobase", "S4Vectors"),
                     ask = FALSE)
install.packages(c("tidyverse", "ggplot2", "ggrepel", "pheatmap", "patchwork", "pander"))
