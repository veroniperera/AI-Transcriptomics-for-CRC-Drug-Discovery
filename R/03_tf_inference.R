source("R/00_setup.R")
suppressPackageStartupMessages({ library(dorothea); library(viper); library(tidyverse) })

deg <- read.csv(file.path(RES, "tables", "DEG_results_full.csv"))
hub <- read.csv(file.path(RES, "tables", "hub_genes_ranked.csv"))

sig_t <- deg %>% filter(!is.na(gene_symbol)) %>% group_by(gene_symbol) %>% slice_max(abs(t), n = 1, with_ties = FALSE) %>% ungroup()
signature <- matrix(sig_t$t, ncol = 1, dimnames = list(sig_t$gene_symbol, "Tumor_vs_Normal"))

data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>% filter(confidence %in% c("A", "B", "C"))
viper_regs <- dorothea2viper_regulons(regulons)

nes <- viper(signature, viper_regs, minsize = 4, eset.filter = FALSE, verbose = FALSE)
tf <- data.frame(TF = rownames(nes), NES = nes[, 1])
tf$p_value <- 2 * pnorm(-abs(tf$NES))
tf$adj_p   <- p.adjust(tf$p_value, method = "BH")
tf$significant <- tf$adj_p < 0.05
tf$is_hub_gene <- tf$TF %in% hub$gene_symbol
tf$module <- hub$module[match(tf$TF, hub$gene_symbol)]
tf <- tf[order(tf$adj_p), ]
write.csv(tf, file.path(RES, "tables", "TableS3_TF_activity.csv"), row.names = FALSE)

tab <- table(Significant = tf$significant, Hub = tf$is_hub_gene)
ft <- fisher.test(tab)
print(tab); cat("Fisher exact p =", ft$p.value, "\n")
overlap <- tf %>% filter(significant, is_hub_gene)
cat("Significant TFs among hub genes:", nrow(overlap), "\n")
write.csv(overlap, file.path(RES, "tables", "TF_hub_overlap.csv"), row.names = FALSE)
writeLines(paste("Fisher exact p =", signif(ft$p.value, 4)), file.path(RES, "tables", "TF_overlap_fisher.txt"))

p <- ggplot(overlap, aes(reorder(TF, NES), NES, fill = module)) +
  geom_col() + coord_flip() + scale_fill_manual(values = c(turquoise = "#00CED1", blue = "#2166AC")) +
  labs(x = NULL, y = "VIPER NES (tumour vs normal)", fill = "Module") + theme_bw(base_size = 12)
ggsave(file.path(RES, "figures", "Fig4A_TF_activity.png"), p, width = 6, height = 5, dpi = 200)
