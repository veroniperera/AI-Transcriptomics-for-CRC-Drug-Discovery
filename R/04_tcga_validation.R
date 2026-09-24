source("R/00_setup.R")
suppressPackageStartupMessages({ library(TCGAbiolinks); library(SummarizedExperiment); library(limma); library(tidyverse) })

query <- GDCquery(project = "TCGA-COAD", data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts")
GDCdownload(query, method = "api", files.per.chunk = 10, directory = file.path(DATA, "raw", "GDCdata"))
se <- GDCprepare(query, directory = file.path(DATA, "raw", "GDCdata"))

info <- as.data.frame(colData(se))
keep_s <- info$sample_type %in% c("Primary Tumor", "Solid Tissue Normal")
se <- se[, keep_s]
cond <- factor(ifelse(colData(se)$sample_type == "Primary Tumor", "Tumor", "Normal"), levels = c("Normal", "Tumor"))
print(table(cond))

counts <- assay(se, "unstranded")
sym <- rowData(se)$gene_name
ord <- order(rowMeans(counts), decreasing = TRUE)
counts <- counts[ord, ]; sym <- sym[ord]
first <- !duplicated(sym) & !is.na(sym)
counts <- counts[first, ]; rownames(counts) <- sym[first]

counts <- counts[rowMeans(counts) > 1, ]
logc <- normalizeBetweenArrays(log2(counts + 1), method = "quantile")
design <- model.matrix(~ cond)
fit <- eBayes(lmFit(logc, design), trend = EBAYES_TREND)
res <- topTable(fit, coef = "condTumor", number = Inf, sort.by = "P")
res$gene_symbol <- rownames(res)
write.csv(res, file.path(RES, "tables", "TCGA_COAD_DEA_full.csv"), row.names = FALSE)

deg <- read.csv(file.path(RES, "tables", "DEG_results_full.csv"))
val <- data.frame(gene_symbol = HUB_TOP10) %>%
  left_join(deg %>% dplyr::select(gene_symbol, logFC_GSE = logFC, adjP_GSE = adj.P.Val), by = "gene_symbol") %>%
  left_join(res %>% dplyr::select(gene_symbol, logFC_TCGA = logFC, adjP_TCGA = adj.P.Val), by = "gene_symbol") %>%
  mutate(same_direction = sign(logFC_GSE) == sign(logFC_TCGA),
         significant_both = adjP_GSE < FDR_CUT & adjP_TCGA < FDR_CUT,
         validated = same_direction & significant_both,
         abs_logFC_TCGA_gt2 = abs(logFC_TCGA) > LOGFC_CUT)
write.csv(val, file.path(RES, "tables", "TableS1_external_validation.csv"), row.names = FALSE)
print(val)
cat("Externally validated hub genes:", sum(val$validated, na.rm = TRUE), "\n")
