library(WGCNA)
library(limma)
library(ggplot2)
library(pheatmap)
library(reshape2)

expr_file <- "expr_file.csv"
header_line <- readLines(expr_file, 1)
sep_char <- ifelse(grepl(",", header_line), ",", "\t")
raw_data <- read.table(expr_file, sep = sep_char, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
rownames(raw_data) <- raw_data[, 1]
expr_data <- raw_data[, -1, drop = FALSE]

dim_names <- list(rownames(expr_data), colnames(expr_data))
expr_num <- matrix(as.numeric(as.matrix(expr_data)), nrow = nrow(expr_data), dimnames = dim_names)
if (max(expr_num, na.rm = TRUE) > 1000) expr_num <- log2(expr_num + 1)
expr_norm <- normalizeBetweenArrays(expr_num)
expr_filtered <- expr_norm[apply(expr_norm, 1, sd, na.rm = TRUE) > 0.5, ]

sample_ids <- colnames(expr_filtered)
groups <- ifelse(grepl("_con$", sample_ids, ignore.case = TRUE), "Control",
                 ifelse(grepl("_tre$", sample_ids, ignore.case = TRUE), "Treatment", "Unknown"))
if (any(groups == "Unknown")) stop("Unrecognized sample names!")

expr_ctrl <- expr_filtered[, which(groups == "Control"), drop = FALSE]
expr_treat <- expr_filtered[, which(groups == "Treatment"), drop = FALSE]
expr_combined <- cbind(expr_ctrl, expr_treat)
expr_matrix <- t(expr_combined)

sample_check <- goodSamplesGenes(expr_matrix, verbose = 3)
if (!sample_check$allOK) {
  expr_matrix <- expr_matrix[sample_check$goodSamples, sample_check$goodGenes]
}

sample_dendro <- hclust(dist(expr_matrix), method = "average")
pdf("Sample_Clustering.pdf", width = 12, height = 9)
par(cex = 0.8, mar = c(0, 4, 2, 0))
plot(sample_dendro, main = "Sample Clustering", sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.2, cex.main = 2)
abline(h = 20000, col = "red", lwd = 2)
dev.off()

cluster_cut <- cutreeStatic(sample_dendro, cutHeight = 20000, minSize = 10)
if (length(unique(cluster_cut)) > 1) {
  expr_matrix <- expr_matrix[cluster_cut == 1, ]
}

trait_data <- data.frame(
  Normal = c(rep(1, sum(groups == "Control")), rep(0, sum(groups == "Treatment"))),
  Disease = c(rep(0, sum(groups == "Control")), rep(1, sum(groups == "Treatment")))
)
rownames(trait_data) <- colnames(expr_combined)

common_samples <- intersect(rownames(expr_matrix), rownames(trait_data))
expr_matrix <- expr_matrix[common_samples, ]
trait_data <- trait_data[common_samples, ]

sample_dendro2 <- hclust(dist(expr_matrix), method = "average")
trait_colors <- numbers2colors(trait_data, signed = FALSE)
pdf("Sample_Heatmap.pdf", width = 12, height = 12)
plotDendroAndColors(sample_dendro2, trait_colors, groupLabels = names(trait_data),
                    main = "Sample Dendrogram and Trait Heatmap", dendroLabels = FALSE)
dev.off()

enableWGCNAThreads()
power_vector <- 1:20
sft_result <- pickSoftThreshold(expr_matrix, powerVector = power_vector, verbose = 5)
optimal_power <- sft_result$powerEstimate

pdf("Scale_Independence_and_Mean_Connectivity.pdf", width = 10, height = 6)
par(mfrow = c(1, 2))
plot(sft_result$fitIndices[, 1], -sign(sft_result$fitIndices[, 3]) * sft_result$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit (signed R^2)",
     type = "n", main = "Scale Independence")
text(sft_result$fitIndices[, 1], -sign(sft_result$fitIndices[, 3]) * sft_result$fitIndices[, 2],
     labels = power_vector, cex = 0.9, col = "blue")
abline(h = 0.90, col = "red", lwd = 2)
plot(sft_result$fitIndices[, 1], sft_result$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = "Mean Connectivity")
text(sft_result$fitIndices[, 1], sft_result$fitIndices[, 5],
     labels = power_vector, cex = 0.9, col = "blue")
dev.off()

adjacency <- adjacency(expr_matrix, power = optimal_power)
tom <- TOMsimilarity(adjacency)
diss_tom <- 1 - tom
gene_dendro <- hclust(as.dist(diss_tom), method = "average")
pdf("Gene_Clustering.pdf", width = 12, height = 9)
plot(gene_dendro, xlab = "", sub = "", main = "Gene Clustering (TOM)", labels = FALSE, hang = 0.04)
dev.off()

dynamic_labels <- cutreeDynamic(dendro = gene_dendro, distM = diss_tom,
                                deepSplit = 2, pamRespectsDendro = FALSE,
                                minClusterSize = 50)
module_colors <- labels2colors(dynamic_labels)
pdf("Dynamic_Tree_Modules.pdf", width = 8, height = 6)
plotDendroAndColors(gene_dendro, module_colors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                    main = "Gene Dendrogram and Module Colors")
dev.off()

me_list <- moduleEigengenes(expr_matrix, colors = module_colors)$eigengenes
me_diss <- 1 - cor(me_list)
me_dendro <- hclust(as.dist(me_diss), method = "average")
merge_threshold <- 0.25
pdf("Module_Clustering.pdf", width = 8, height = 6)
plot(me_dendro, main = "Clustering of Module Eigengenes", xlab = "", sub = "")
abline(h = merge_threshold, col = "red", lwd = 2)
dev.off()

merge_result <- mergeCloseModules(expr_matrix, module_colors, cutHeight = merge_threshold, verbose = 3)
merged_colors <- merge_result$colors
merged_mes <- merge_result$newMEs

pdf("Merged_Modules_Comparison.pdf", width = 10, height = 6)
plotDendroAndColors(gene_dendro, cbind(module_colors, merged_colors),
                    c("Original", "Merged"), dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, main = "Module Comparison")
dev.off()

module_colors <- merged_colors
final_mes <- merged_mes

mod_trait_cor <- cor(final_mes, trait_data, use = "p")
mod_trait_p <- corPvalueStudent(mod_trait_cor, nrow(expr_matrix))

cor_melt <- melt(mod_trait_cor)
pval_melt <- melt(mod_trait_p)
heatmap_df <- merge(cor_melt, pval_melt, by = c("Var1", "Var2"))
colnames(heatmap_df) <- c("Module", "Trait", "Correlation", "Pvalue")

pdf("Module_Trait_Heatmap.pdf", width = 5, height = 4)
ggplot(heatmap_df, aes(x = Trait, y = Module, fill = Correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#4DBBD5", mid = "#F7F7F7", high = "#DC0000",
                       midpoint = 0, limits = c(-1, 1)) +
  geom_text(aes(label = sprintf("%.2f\n(%.1e)", Correlation, Pvalue)), size = 4.6) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16)) +
  ggtitle("Module-Trait Correlation")
dev.off()

gene_mm <- as.data.frame(cor(expr_matrix, final_mes, use = "p"))
gene_mm_p <- as.data.frame(corPvalueStudent(as.matrix(gene_mm), nrow(expr_matrix)))
gene_gs <- as.data.frame(cor(expr_matrix, trait_data, use = "p"))
gene_gs_p <- as.data.frame(corPvalueStudent(as.matrix(gene_gs), nrow(expr_matrix)))

gene_names <- colnames(expr_matrix)
gene_info <- data.frame(Gene = gene_names, Module = module_colors)
for (mod in colnames(final_mes)) {
  gene_info[, paste0("MM_", mod)] <- gene_mm[, mod]
  gene_info[, paste0("p.MM_", mod)] <- gene_mm_p[, mod]
}
for (trait in colnames(trait_data)) {
  gene_info[, paste0("GS_", trait)] <- gene_gs[, trait]
  gene_info[, paste0("p.GS_", trait)] <- gene_gs_p[, trait]
}
gene_info <- gene_info[order(gene_info$Module), ]
write.table(gene_info, file = "GeneInfo_Modules.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

unique_modules <- unique(module_colors)
for (mod in unique_modules) {
  module_genes <- gene_names[module_colors == mod]
  write.table(module_genes, file = paste0("ModuleGenes_", mod, ".txt"),
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

module_sizes <- table(module_colors)
module_size_df <- as.data.frame(module_sizes)
colnames(module_size_df) <- c("Module", "GeneCount")
pdf("Module_Gene_Counts.pdf", width = 8, height = 6)
ggplot(module_size_df, aes(x = Module, y = GeneCount, fill = Module)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
                               "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85")) +
  theme_minimal() + ggtitle("Gene Counts per Module") + xlab("Module") + ylab("Gene Count") +
  theme(legend.position = "none")
dev.off()
