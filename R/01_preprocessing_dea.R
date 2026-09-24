source("R/00_setup.R")
suppressPackageStartupMessages({
  library(GEOquery); library(limma); library(tidyverse); library(ggrepel); library(org.Hs.eg.db)
})

gse <- getGEO("GSE156451", GSEMatrix = TRUE, getGPL = FALSE)
meta_all <- pData(gse[[1]])
meta <- meta_all %>%
  dplyr::select(geo_accession, characteristics_ch1, `tissue:ch1`) %>%
  dplyr::rename(tissue = `tissue:ch1`)

getGEOSuppFiles("GSE156451", baseDir = file.path(DATA, "raw"))
count_file <- list.files(file.path(DATA, "raw", "GSE156451"), pattern = "count|Count|CPM|TMM", full.names = TRUE)[1]
expr_all <- read.table(count_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

tumor  <- meta %>% filter(tissue == "Tumor")         %>% head(25) %>% pull(geo_accession)
normal <- meta %>% filter(tissue == "Native tissue") %>% head(25) %>% pull(geo_accession)
sel    <- c(tumor, normal)
meta   <- meta %>% filter(geo_accession %in% sel)
meta   <- meta[match(sel, meta$geo_accession), ]
counts <- expr_all[, sel]
stopifnot(sum(is.na(counts)) == 0)

keep <- rowMeans(counts) > 1
counts_f <- counts[keep, ]
cat("Genes retained:", nrow(counts_f), "\n")

group  <- factor(meta$tissue, levels = c("Native tissue", "Tumor"))
design <- model.matrix(~ group)
colnames(design) <- c("Intercept", "Tumor_vs_Native")

log_raw <- log2(as.matrix(counts_f) + 1)
log_expr <- if (QUANTILE_NORM) normalizeBetweenArrays(log_raw, method = "quantile") else log_raw
dimnames(log_expr) <- dimnames(log_raw)

fit <- eBayes(lmFit(log_expr, design), trend = EBAYES_TREND)
res <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
res$entrez_id <- rownames(res)
sym <- mapIds(org.Hs.eg.db, keys = rownames(res), column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
nm  <- mapIds(org.Hs.eg.db, keys = rownames(res), column = "GENENAME", keytype = "ENTREZID", multiVals = "first")
res$gene_symbol <- ifelse(is.na(sym), rownames(res), sym)
res$gene_name   <- ifelse(is.na(nm), "", nm)
res$status <- "Not significant"
res$status[res$logFC >  LOGFC_CUT & res$adj.P.Val < FDR_CUT] <- "Upregulated"
res$status[res$logFC < -LOGFC_CUT & res$adj.P.Val < FDR_CUT] <- "Downregulated"
res <- res[, c("entrez_id", "gene_symbol", "gene_name", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "status")]

up_n <- sum(res$status == "Upregulated"); down_n <- sum(res$status == "Downregulated")
cat("DEGs:", up_n + down_n, "| up:", up_n, "| down:", down_n, "\n")

write.csv(res, file.path(RES, "tables", "DEG_results_full.csv"), row.names = FALSE)
write.csv(res[res$status != "Not significant", ], file.path(RES, "tables", "DEG_results_significant.csv"), row.names = FALSE)
write.csv(meta, file.path(DATA, "processed", "metadata.csv"), row.names = FALSE)
saveRDS(log_expr, file.path(RES, "rds", "log_expr.rds"))
saveRDS(meta,     file.path(RES, "rds", "meta.rds"))

cols <- c(Tumor = "#E05C5C", `Native tissue` = "#5B9BD5")
idx <- c(which(group == "Tumor")[1:5], which(group == "Native tissue")[1:5])
png(file.path(RES, "figures", "Fig1_normalization.png"), width = 1800, height = 850, res = 130)
par(mfrow = c(1, 2), mar = c(9, 5, 4, 2))
boxplot(log_raw[, idx], col = cols[as.character(group[idx])], las = 2, cex.axis = 0.7,
        ylab = "log2(counts + 1)", main = "Before normalization", outline = FALSE)
boxplot(log_expr[, idx], col = cols[as.character(group[idx])], las = 2, cex.axis = 0.7,
        ylab = "Normalized log2 expression", main = "After normalization", outline = FALSE)
dev.off()

pca <- prcomp(t(log_expr), scale. = FALSE)
pv  <- round(summary(pca)$importance[2, 1:2] * 100, 1)
p_pca <- ggplot(data.frame(pca$x[, 1:2], condition = group), aes(PC1, PC2, color = condition)) +
  geom_point(size = 3.5, alpha = 0.9) + scale_color_manual(values = cols) +
  labs(x = paste0("PC1 (", pv[1], "%)"), y = paste0("PC2 (", pv[2], "%)"), color = "Condition") + theme_bw(base_size = 13)
ggsave(file.path(RES, "figures", "Fig_PCA.png"), p_pca, width = 7, height = 6, dpi = 150)

lab <- res %>% filter(status != "Not significant") %>% arrange(adj.P.Val) %>% head(15)
p_vol <- ggplot(res, aes(logFC, -log10(adj.P.Val), color = status)) +
  geom_point(size = 0.9, alpha = 0.7) +
  scale_color_manual(values = c(Upregulated = "#E05C5C", Downregulated = "#5B9BD5", `Not significant` = "grey75")) +
  geom_vline(xintercept = c(-LOGFC_CUT, LOGFC_CUT), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(FDR_CUT), linetype = "dashed", color = "grey40") +
  geom_text_repel(data = lab, aes(label = gene_symbol), size = 2.8, max.overlaps = 20, show.legend = FALSE) +
  labs(x = "log2 fold change", y = "-log10 adjusted p-value", color = "DEG status") + theme_bw(base_size = 13)
ggsave(file.path(RES, "figures", "Fig_Volcano.png"), p_vol, width = 8, height = 7, dpi = 150)
