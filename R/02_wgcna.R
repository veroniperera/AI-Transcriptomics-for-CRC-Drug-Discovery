source("R/00_setup.R")
suppressPackageStartupMessages({ library(WGCNA); library(tidyverse); library(org.Hs.eg.db) })
allowWGCNAThreads()

log_expr <- readRDS(file.path(RES, "rds", "log_expr.rds"))
meta     <- readRDS(file.path(RES, "rds", "meta.rds"))
deg      <- read.csv(file.path(RES, "tables", "DEG_results_full.csv"), colClasses = c(entrez_id = "character"))

if (!is.null(WGCNA_TOP_MAD)) {
  mads <- apply(log_expr, 1, mad)
  log_expr <- log_expr[names(sort(mads, decreasing = TRUE))[seq_len(WGCNA_TOP_MAD)], , drop = FALSE]
}
datExpr <- t(log_expr)
gsg <- goodSamplesGenes(datExpr, verbose = 0)
datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
cat("WGCNA input:", dim(datExpr), "\n")

traits <- data.frame(TumorStatus = as.numeric(meta$tissue == "Tumor"),
                     NormalStatus = as.numeric(meta$tissue != "Tumor"),
                     row.names = rownames(datExpr))

sft <- pickSoftThreshold(datExpr, powerVector = 1:20, networkType = "unsigned", verbose = 0)
write.csv(sft$fitIndices, file.path(RES, "tables", "soft_threshold_fit.csv"), row.names = FALSE)

net <- blockwiseModules(datExpr, power = SOFT_POWER, networkType = "unsigned", TOMType = "unsigned",
                        minModuleSize = 30, deepSplit = 2, mergeCutHeight = 0.25,
                        maxBlockSize = ncol(datExpr) + 1, numericLabels = FALSE,
                        pamRespectsDendro = FALSE, saveTOMs = FALSE, randomSeed = SEED, verbose = 0)
MEs <- orderMEs(net$MEs)
cat("Modules:", length(unique(net$colors)), "\n")

ME_cor <- cor(MEs, traits, use = "p")
ME_p   <- corPvalueStudent(ME_cor, nrow(datExpr))
mt <- data.frame(module = rownames(ME_cor), r_tumour = ME_cor[, "TumorStatus"], p_tumour = ME_p[, "TumorStatus"])
mt <- mt[order(-abs(mt$r_tumour)), ]
write.csv(mt, file.path(RES, "tables", "module_trait_correlation.csv"), row.names = FALSE)

kept <- mt[abs(mt$r_tumour) > ME_TRAIT_R & mt$p_tumour < 0.05 & mt$module != "MEgrey", ]
top2 <- sub("^ME", "", head(kept$module, 2))
cat("Modules carried forward:", top2, "\n")

GS <- cor(datExpr, traits$TumorStatus, use = "p")[, 1]
hub <- do.call(rbind, lapply(top2, function(mod) {
  genes <- names(net$colors)[net$colors == mod]
  MM <- cor(datExpr[, genes], MEs[, paste0("ME", mod)], use = "p")[, 1]
  data.frame(entrez_id = genes, module = mod, MM = MM, GS = GS[genes])
}))
hub$hub <- abs(hub$MM) > MM_CUT & abs(hub$GS) > GS_CUT
hub <- hub[hub$hub, ]
hub$gene_symbol <- deg$gene_symbol[match(hub$entrez_id, deg$entrez_id)]
hub$logFC       <- deg$logFC[match(hub$entrez_id, deg$entrez_id)]
hub$adj.P.Val   <- deg$adj.P.Val[match(hub$entrez_id, deg$entrez_id)]
hub$is_DEG      <- deg$status[match(hub$entrez_id, deg$entrez_id)] != "Not significant"
hub$rank_score  <- rank(-abs(hub$MM)) + rank(-abs(hub$GS))
hub <- hub[order(hub$rank_score), ]
write.csv(hub, file.path(RES, "tables", "hub_genes_ranked.csv"), row.names = FALSE)
print(table(hub$module)); cat("Hub genes that are also DEGs:", sum(hub$is_DEG), "\n")

for (mod in top2) {
  h <- hub[hub$module == mod, ]
  cat(mod, "| n =", nrow(h), "| MM<0:", sum(h$MM < 0), "| GS<0:", sum(h$GS < 0), "\n")
}

saveRDS(list(net = net, MEs = MEs, traits = traits, sft = sft, datExpr = datExpr), file.path(RES, "rds", "wgcna.rds"))

png(file.path(RES, "figures", "Fig3A_soft_threshold.png"), width = 2000, height = 900, res = 200)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
r2 <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
plot(sft$fitIndices[, 1], r2, type = "n", xlab = "Soft threshold (power)", ylab = "Scale-free topology R2", ylim = c(0, 1))
text(sft$fitIndices[, 1], r2, labels = sft$fitIndices[, 1], col = "#E05C5C", cex = 0.8)
abline(h = 0.85, col = "grey40", lty = 2); abline(v = SOFT_POWER, col = "#E05C5C", lty = 2)
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], type = "n", xlab = "Soft threshold (power)", ylab = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = sft$fitIndices[, 1], col = "#5B9BD5", cex = 0.8)
abline(v = SOFT_POWER, col = "#E05C5C", lty = 2)
dev.off()

png(file.path(RES, "figures", "Fig3B_dendrogram.png"), width = 2400, height = 1100, res = 200)
plotDendroAndColors(net$dendrograms[[1]], net$colors[net$blockGenes[[1]]], "Module colours",
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "")
dev.off()

textMatrix <- matrix(paste0(signif(ME_cor, 2), "\n(p=", signif(ME_p, 1), ")"), nrow = nrow(ME_cor))
png(file.path(RES, "figures", "Fig3C_module_trait.png"), width = 1500, height = 2000, res = 200)
par(mar = c(4, 8, 1, 2))
labeledHeatmap(Matrix = ME_cor, xLabels = c("Tumour status", "Normal status"), yLabels = rownames(ME_cor),
               ySymbols = rownames(ME_cor), colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.55, main = "")
dev.off()

pdf(file.path(RES, "figures", "Fig3D_MM_vs_GS.pdf"), width = 9, height = 4.5)
par(mfrow = c(1, 2))
for (mod in top2) {
  genes <- names(net$colors)[net$colors == mod]
  MM <- abs(cor(datExpr[, genes], MEs[, paste0("ME", mod)], use = "p")[, 1])
  g  <- abs(GS[genes])
  verboseScatterplot(MM, g, xlab = paste("Module membership -", mod), ylab = "Gene significance - tumour status",
                     col = ifelse(MM > MM_CUT & g > GS_CUT, mod, "grey80"), abline = TRUE)
  abline(v = MM_CUT, h = GS_CUT, lty = 2, col = "grey50")
}
dev.off()
