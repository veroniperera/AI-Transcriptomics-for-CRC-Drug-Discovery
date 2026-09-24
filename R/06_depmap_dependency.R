source("R/00_setup.R")
suppressPackageStartupMessages(library(tidyverse))

effect_file <- file.path(DATA, "raw", "CRISPRGeneEffect.csv")
stopifnot(file.exists(effect_file))
DEP_THRESHOLD <- -0.5

eff <- read.csv(effect_file, check.names = FALSE, row.names = 1)
colnames(eff) <- sub(" \\(.*\\)$", "", colnames(eff))

summ <- do.call(rbind, lapply(c("SHMT2", "RNF19A"), function(g) {
  x <- eff[[g]]
  data.frame(gene = g, n_lines = sum(!is.na(x)), n_dependent = sum(x <= DEP_THRESHOLD, na.rm = TRUE),
             mean_effect = mean(x, na.rm = TRUE), median_effect = median(x, na.rm = TRUE))
}))
print(summ)
write.csv(summ, file.path(RES, "tables", "Table4B_depmap_dependency.csv"), row.names = FALSE)
