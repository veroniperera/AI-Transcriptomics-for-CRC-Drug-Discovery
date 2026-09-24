source("R/00_setup.R")
suppressPackageStartupMessages({
  library(PharmacoGx); library(Biobase); library(SummarizedExperiment); library(tidyverse)
})

rds <- file.path(DATA, "raw", "GDSC.rds")
if (file.exists(rds)) GDSC <- readRDS(rds) else { GDSC <- downloadPSet("GDSC_2020(v2-8.2)"); saveRDS(GDSC, rds) }

cell_info <- cellInfo(GDSC)
crc_lines <- rownames(cell_info)[grep("Bowel", cell_info$tissueid, ignore.case = TRUE)]
cat("CRC lines:", length(crc_lines), "\n")

targets  <- c("SHMT2", "RNF19A")
feat_ids <- fNames(GDSC, "rna")[which(featureInfo(GDSC, "rna")$Symbol %in% targets)]
feat_sym <- featureInfo(GDSC, "rna")[feat_ids, "Symbol"]
stopifnot(length(feat_ids) == 2)

expr <- assay(summarizeMolecularProfiles(GDSC, mDataType = "rna", features = feat_ids, verbose = FALSE), 1)
rownames(expr) <- feat_sym
expr_crc <- expr[, colnames(expr) %in% crc_lines]
write.csv(expr_crc, file.path(RES, "tables", "GDSC_CRC_expression_SHMT2_RNF19A.csv"))
cat("Mean log2 expression in CRC lines:\n"); print(rowMeans(expr_crc, na.rm = TRUE))

ic50 <- summarizeSensitivityProfiles(GDSC, sensitivity.measure = "ic50_recomputed", summary.stat = "median", verbose = FALSE)
ic50_crc <- ic50[, colnames(ic50) %in% crc_lines]
cat("Drugs:", nrow(ic50_crc), "| CRC lines:", ncol(ic50_crc), "\n")

sig <- drugSensitivitySig(object = GDSC, mDataType = "rna", drugs = drugNames(GDSC), features = feat_ids,
                          sensitivity.measure = "ic50_recomputed", molecular.summary.stat = "median",
                          sensitivity.summary.stat = "median", verbose = FALSE)
saveRDS(sig, file.path(RES, "rds", "PharmacoGx_sig_results.rds"))
fdr <- sig[, , "fdr"]; est <- sig[, , "estimate"]; tst <- sig[, , "tstat"]
rownames(fdr) <- rownames(est) <- rownames(tst) <- feat_sym

tab <- do.call(rbind, lapply(targets, function(g) {
  data.frame(gene = g, drug = colnames(est), estimate = est[g, ], tstat = tst[g, ], fdr = fdr[g, ])
}))
tab <- tab[order(tab$gene, tab$fdr), ]
write.csv(tab, file.path(RES, "tables", "TableS4_gene_drug_associations_signature.csv"), row.names = FALSE)
print(tab %>% filter(fdr < FDR_CUT) %>% group_by(gene) %>% summarise(n_sig = n(), top_drug = drug[1], top_fdr = fdr[1]))

cl <- intersect(colnames(expr_crc), colnames(ic50_crc))
crc_assoc <- do.call(rbind, lapply(targets, function(g) {
  do.call(rbind, lapply(rownames(ic50_crc), function(d) {
    x <- as.numeric(expr_crc[g, cl]); y <- log2(as.numeric(ic50_crc[d, cl]))
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 10) return(NULL)
    ct <- cor.test(x[ok], y[ok], method = "pearson")
    data.frame(gene = g, drug = d, n = sum(ok), r = unname(ct$estimate), p = ct$p.value)
  }))
}))
crc_assoc$fdr <- ave(crc_assoc$p, crc_assoc$gene, FUN = function(p) p.adjust(p, method = "BH"))
write.csv(crc_assoc[order(crc_assoc$gene, crc_assoc$fdr), ],
          file.path(RES, "tables", "TableS4b_CRC_only_pearson.csv"), row.names = FALSE)

long <- as.data.frame(t(expr_crc)) %>% rownames_to_column("CellLine") %>%
  pivot_longer(all_of(targets), names_to = "Gene", values_to = "Expression")
p1 <- ggplot(long, aes(Gene, Expression, fill = Gene)) + geom_boxplot(outlier.shape = 21, width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.5) +
  scale_fill_manual(values = c(RNF19A = "#2166AC", SHMT2 = "#B2182B")) +
  labs(x = NULL, y = "log2 expression") + theme_bw(base_size = 12) + theme(legend.position = "none")
ggsave(file.path(RES, "figures", "Fig4B_expression_boxplot.png"), p1, width = 5, height = 5, dpi = 300)

waterfall <- function(g, col) {
  d <- tab %>% filter(gene == g, !is.na(estimate)) %>% arrange(estimate) %>%
    mutate(drug = factor(drug, levels = drug), sig = fdr < FDR_CUT)
  ggplot(d, aes(drug, estimate, fill = sig)) + geom_col() +
    scale_fill_manual(values = c(`TRUE` = col, `FALSE` = "grey70"), labels = c("FDR >= 0.05", "FDR < 0.05")) +
    geom_hline(yintercept = 0) + labs(x = NULL, y = "Association estimate", fill = NULL, title = g) +
    theme_bw(base_size = 11) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
}
ggsave(file.path(RES, "figures", "Fig4C_waterfall_RNF19A.png"), waterfall("RNF19A", "#B2182B"), width = 8, height = 5, dpi = 300)
ggsave(file.path(RES, "figures", "Fig4D_waterfall_SHMT2.png"), waterfall("SHMT2", "#2166AC"), width = 8, height = 5, dpi = 300)

scatter <- function(g, d, col) {
  df <- data.frame(expr = as.numeric(expr_crc[g, cl]), ic50 = as.numeric(ic50_crc[d, cl]))
  ggplot(na.omit(df), aes(expr, ic50)) + geom_point(color = col, size = 2.5, alpha = 0.8) +
    geom_smooth(method = "lm", color = "grey30") +
    labs(x = paste(g, "log2 expression"), y = paste(d, "IC50 (recomputed)")) + theme_bw(base_size = 12)
}
ggsave(file.path(RES, "figures", "Fig4E_RNF19A_Sorafenib.png"), scatter("RNF19A", "Sorafenib", "#2166AC"), width = 5, height = 5, dpi = 300)
ggsave(file.path(RES, "figures", "Fig4F_SHMT2_Sorafenib.png"), scatter("SHMT2", "Sorafenib", "#B2182B"), width = 5, height = 5, dpi = 300)
