library(reshape2)
library(WGCNA)
library(limma)
library(pheatmap)
library(ggplot2)
library(grid)
library(RColorBrewer)

set.seed(1234)

workDir <- "C:"
expFilePath <- "merged.csv"

if (!dir.exists(workDir)) stop("Working directory does not exist!")
setwd(workDir)
if (!file.exists(expFilePath)) stop("Expression matrix file not found!")

output_dir <- file.path(workDir, "WGCNA_output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
cat("Output directory:", output_dir, "\n")

nature_gradient <- c("#053061", "#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B", "#67001F")
nature_colors   <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
                     "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85")

mod_colors_palette <- c(
  "brown"="#8B4513", "blue"="#4DBBD5", "yellow"="#FFC107", "turquoise"="#00A087",
  "grey"="#A9A9A9", "black"="#333333", "magenta"="#E64B35",
  "green"="#00A087", "red"="#DC0000", "pink"="#F39B7F",
  "purple"="#8491B4", "salmon"="#FA8072", "orange"="#FFA500",
  "cyan"="#91D1C2", "lightyellow"="#FFFACD",
  "darkorange"="#FF8C00", "royalblue"="#3C5488", "darkgrey"="#7E6148"
)

tmp_head <- readLines(expFilePath, 1)
sep <- ifelse(grepl(",", tmp_head), ",", "\t")
rawData <- read.table(expFilePath, sep = sep, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
rownames(rawData) <- rawData[, 1]
exprData <- rawData[, -1]

dimNames <- list(rownames(exprData), colnames(exprData))
numericData <- matrix(as.numeric(as.matrix(exprData)), nrow = nrow(exprData), dimnames = dimNames)
if (max(numericData, na.rm = TRUE) > 1000) numericData <- log2(numericData + 1)
normalizedData <- normalizeBetweenArrays(numericData)
filteredData <- normalizedData[apply(normalizedData, 1, sd, na.rm = TRUE) > 0.5, ]

sample_names <- colnames(filteredData)
group_info <- ifelse(grepl("_con$", sample_names, ignore.case = TRUE), "Control",
                     ifelse(grepl("_tre$", sample_names, ignore.case = TRUE), "Treatment", "Unknown"))
if (any(group_info == "Unknown")) stop("Unrecognized sample names!")

dataControl   <- filteredData[, which(group_info == "Control"), drop=FALSE]
dataTreatment <- filteredData[, which(group_info == "Treatment"), drop=FALSE]
combinedData  <- cbind(dataControl, dataTreatment)
exprMatrix    <- t(combinedData)

sampleCheck <- goodSamplesGenes(exprMatrix, verbose = 3)
if (!sampleCheck$allOK) {
  exprMatrix <- exprMatrix[sampleCheck$goodSamples, sampleCheck$goodGenes]
}

sampleDendro <- hclust(dist(exprMatrix), method = "average")
pdf(file.path(output_dir, "Sample_Clustering.pdf"), width = 12, height = 9)
par(cex = 0.8, mar = c(0, 4, 2, 0))
plot(sampleDendro, main = "Sample Clustering", sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.2, cex.main = 2)
abline(h = 20000, col = "red", lwd = 2)
dev.off()

clusterCut <- cutreeStatic(sampleDendro, cutHeight = 20000, minSize = 10)
if (length(unique(clusterCut)) > 1) {
  exprMatrix <- exprMatrix[clusterCut == 1, ]
}

clinicalData <- data.frame(
  Normal = c(rep(1, sum(group_info == "Control")), rep(0, sum(group_info == "Treatment"))),
  Disease = c(rep(0, sum(group_info == "Control")), rep(1, sum(group_info == "Treatment")))
)
rownames(clinicalData) <- colnames(combinedData)

commonSamples <- intersect(rownames(exprMatrix), rownames(clinicalData))
exprMatrix <- exprMatrix[commonSamples, ]
clinicalData <- clinicalData[commonSamples, ]

sampleDendro2 <- hclust(dist(exprMatrix), method = "average")
traitColors <- numbers2colors(clinicalData, signed = FALSE)
pdf(file.path(output_dir, "Sample_Heatmap.pdf"), width = 12, height = 12)
plotDendroAndColors(sampleDendro2, traitColors, groupLabels = names(clinicalData),
                    main = "Sample Dendrogram and Trait Heatmap", dendroLabels = FALSE)
dev.off()

closeAllConnections()
enableWGCNAThreads()

powerVector <- 1:20
sftResult <- pickSoftThreshold(exprMatrix, powerVector = powerVector, verbose = 5)
optimalPower <- sftResult$powerEstimate
cat("Optimal soft threshold power:", optimalPower, "\n")

pdf(file.path(output_dir, "Scale_Independence_and_Mean_Connectivity.pdf"), width = 10, height = 6)
par(mfrow = c(1, 2))
plot(sftResult$fitIndices[, 1], -sign(sftResult$fitIndices[, 3]) * sftResult$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit (signed R^2)",
     type = "n", main = "Scale Independence")
text(sftResult$fitIndices[, 1], -sign(sftResult$fitIndices[, 3]) * sftResult$fitIndices[, 2],
     labels = powerVector, cex = 0.9, col = "blue")
abline(h = 0.90, col = "red", lwd = 2)
plot(sftResult$fitIndices[, 1], sftResult$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = "Mean Connectivity")
text(sftResult$fitIndices[, 1], sftResult$fitIndices[, 5],
     labels = powerVector, cex = 0.9, col = "blue")
dev.off()

adjacencyMatrix <- adjacency(exprMatrix, power = optimalPower)
TOMMatrix <- TOMsimilarity(adjacencyMatrix)
dissTOM <- 1 - TOMMatrix

geneDendro <- hclust(as.dist(dissTOM), method = "average")
pdf(file.path(output_dir, "Gene_Clustering.pdf"), width = 12, height = 9)
plot(geneDendro, xlab = "", sub = "", main = "Gene Clustering (TOM)", labels = FALSE, hang = 0.04)
dev.off()

dynamicModuleLabels <- cutreeDynamic(dendro = geneDendro, distM = dissTOM,
                                     deepSplit = 2, pamRespectsDendro = FALSE,
                                     minClusterSize = 50)
moduleColors <- labels2colors(dynamicModuleLabels)
pdf(file.path(output_dir, "Dynamic_Tree_Modules.pdf"), width = 8, height = 6)
plotDendroAndColors(geneDendro, moduleColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                    main = "Gene Dendrogram and Module Colors")
dev.off()

MEs <- moduleEigengenes(exprMatrix, colors = moduleColors)$eigengenes
moduleDiss <- 1 - cor(MEs)
moduleEigDendro <- hclust(as.dist(moduleDiss), method = "average")
mergeThreshold <- 0.25
pdf(file.path(output_dir, "Module_Clustering.pdf"), width = 8, height = 6)
plot(moduleEigDendro, main = "Clustering of Module Eigengenes", xlab = "", sub = "")
abline(h = mergeThreshold, col = "red", lwd = 2)
dev.off()

mergeResult <- mergeCloseModules(exprMatrix, moduleColors, cutHeight = mergeThreshold, verbose = 3)
mergedModuleColors <- mergeResult$colors
mergedMEs <- mergeResult$newMEs

pdf(file.path(output_dir, "Merged_Modules_Comparison.pdf"), width = 10, height = 6)
plotDendroAndColors(geneDendro, cbind(moduleColors, mergedModuleColors),
                    c("Original", "Merged"), dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, main = "Module Comparison")
dev.off()

moduleColors <- mergedModuleColors
moduleEigengenes <- mergedMEs

moduleTraitCor <- cor(moduleEigengenes, clinicalData, use = "p")
moduleTraitPvalues <- corPvalueStudent(moduleTraitCor, nrow(exprMatrix))

corDF <- melt(moduleTraitCor)
pvalDF <- melt(moduleTraitPvalues)
heatmapDF <- merge(corDF, pvalDF, by = c("Var1", "Var2"))
colnames(heatmapDF) <- c("Module", "Trait", "Correlation", "Pvalue")

pdf(file.path(output_dir, "Module_Trait_Heatmap.pdf"), width = 5, height = 4)
ggplot(heatmapDF, aes(x = Trait, y = Module, fill = Correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#4DBBD5", mid = "#F7F7F7", high = "#DC0000",
                       midpoint = 0, limits = c(-1, 1)) +
  geom_text(aes(label = sprintf("%.2f\n(%.1e)", Correlation, Pvalue)), size = 4.6) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16)) +
  ggtitle("Module-Trait Correlation")
dev.off()

geneModuleMembership <- as.data.frame(cor(exprMatrix, moduleEigengenes, use = "p"))
geneMMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(exprMatrix)))
geneTraitSignificance <- as.data.frame(cor(exprMatrix, clinicalData, use = "p"))
geneGSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nrow(exprMatrix)))

geneNames <- colnames(exprMatrix)
geneInfo <- data.frame(Gene = geneNames, Module = moduleColors)
for (mod in colnames(moduleEigengenes)) {
  geneInfo[, paste0("MM_", mod)] <- geneModuleMembership[, mod]
  geneInfo[, paste0("p.MM_", mod)] <- geneMMPvalue[, mod]
}
for (trait in colnames(clinicalData)) {
  geneInfo[, paste0("GS_", trait)] <- geneTraitSignificance[, trait]
  geneInfo[, paste0("p.GS_", trait)] <- geneGSPvalue[, trait]
}
geneInfo <- geneInfo[order(geneInfo$Module), ]
write.table(geneInfo, file = file.path(output_dir, "GeneInfo_Modules.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

uniqueModules <- unique(moduleColors)
for (mod in uniqueModules) {
  moduleGenes <- geneNames[moduleColors == mod]
  write.table(moduleGenes, file = file.path(output_dir, paste0("ModuleGenes_", mod, ".txt")),
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

moduleSizes <- table(moduleColors)
moduleSizeDF <- as.data.frame(moduleSizes)
colnames(moduleSizeDF) <- c("Module", "GeneCount")
pdf(file.path(output_dir, "Module_Gene_Counts.pdf"), width = 8, height = 6)
ggplot(moduleSizeDF, aes(x = Module, y = GeneCount, fill = Module)) +
  geom_bar(stat = "identity") + scale_fill_manual(values = nature_colors) +
  theme_minimal() + ggtitle("Gene Counts per Module") + xlab("Module") + ylab("Gene Count") +
  theme(legend.position = "none")
dev.off()

targetTrait <- "Disease"
for (mod in unique(moduleColors)) {
  pdf(file = file.path(output_dir, paste0("MM_vs_GS_", mod, ".pdf")), width = 6, height = 6)
  inModule <- (moduleColors == mod)
  mmCol <- paste0("ME", mod)
  MM <- as.numeric(geneModuleMembership[inModule, mmCol])
  GS <- as.numeric(geneTraitSignificance[inModule, targetTrait])
  ct <- cor.test(MM, GS)
  plot(MM, GS,
       xlab = paste("Module Membership in", mod),
       ylab = paste("Gene Significance for", targetTrait),
       main = paste0(mod, " cor=", signif(ct$estimate, 3), ", p=", format(ct$p.value, scientific = TRUE, digits = 2)),
       pch = 21, bg = adjustcolor(mod, alpha.f = 0.6), col = "black", cex = 1.5)
  abline(lm(GS ~ MM), col = "blue", lwd = 2, lty = 2)
  abline(v = 0.8, h = 0.2, col = "orange", lty = 3, lwd = 1.5)
  dev.off()
}

df <- data.frame(
  Module = moduleColors,
  MM = sapply(seq_along(moduleColors), function(i) geneModuleMembership[i, paste0("ME", moduleColors[i])]),
  GS = geneTraitSignificance[, "Disease"]
)
pdf(file.path(output_dir, "Boxplot_MM.pdf"), width = 7, height = 6)
ggplot(df, aes(x = Module, y = MM, fill = Module)) +
  geom_boxplot(alpha = 0.7) + scale_fill_manual(values = nature_colors) +
  theme_minimal() + labs(title = "Module Membership per Module", y = "Module Membership") +
  theme(legend.position = "none")
dev.off()

pdf(file.path(output_dir, "ViolinPlot_GS.pdf"), width = 7, height = 6)
ggplot(df, aes(x = Module, y = GS, fill = Module)) +
  geom_violin(alpha = 0.7) + scale_fill_manual(values = nature_colors) +
  theme_minimal() + labs(title = "Gene Significance per Module", y = "Gene Significance") +
  theme(legend.position = "none")
dev.off()

moduleCor <- cor(moduleEigengenes, use = "p")
moduleCor_melt <- melt(moduleCor)
pdf(file.path(output_dir, "Module_Module_Correlation.pdf"), width = 7, height = 6)
ggplot(moduleCor_melt, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") + scale_fill_gradientn(colors = nature_gradient, limits = c(-1, 1)) +
  geom_text(aes(label = sprintf("%.2f", value)), size = 3) +
  theme_minimal() + ggtitle("Module-Module Correlation") +
  xlab("Module Eigengenes") + ylab("Module Eigengenes") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
dev.off()

targetModule <- "green"
inModule <- (moduleColors == targetModule)
if (sum(inModule) > 1) {
  modGenes <- geneNames[inModule]
  modTOM <- TOMMatrix[inModule, inModule]
  geneDendroModule <- hclust(as.dist(1 - modTOM), method = "average")
  pdf(file.path(output_dir, paste0("TOMplot_", targetModule, ".pdf")), width = 10, height = 10)
  TOMplot(modTOM, geneDendroModule, main = paste("TOM Plot for", targetModule, "Module"))
  dev.off()
}

for (targetMod in unique(moduleColors)) {
  mod_genes <- geneNames[moduleColors == targetMod]
  if (length(mod_genes) < 2) next
  mod_color <- ifelse(targetMod %in% names(mod_colors_palette), mod_colors_palette[targetMod], "#888888")
  me_col <- paste0("ME", targetMod)
  if (!(me_col %in% colnames(MEs))) next
  sample_list <- rownames(exprMatrix)
  mod_expr <- exprMatrix[sample_list, mod_genes, drop = FALSE]
  eigengene <- MEs[sample_list, me_col]
  bar_df <- data.frame(Sample = factor(sample_list, levels = sample_list), Eigengene = eigengene)
  
  heatmap_grob <- grid.grabExpr({
    pheatmap(t(mod_expr), cluster_cols = FALSE, cluster_rows = TRUE,
             show_rownames = FALSE, show_colnames = FALSE, scale = "row",
             color = colorRampPalette(nature_gradient)(100),
             border_color = NA, legend = FALSE, fontsize = 10)
  })
  
  pdf_file <- file.path(output_dir, paste0("Module_", targetMod, "_heatmap_eigengene.pdf"))
  pdf(pdf_file, width = 9, height = 7)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(c(0.58, 0.38), "npc"))))
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid.draw(heatmap_grob)
  grid.text(targetMod, x = 0.5, y = unit(0.98, "npc"),
            gp = gpar(fontface = "bold", fontsize = 20, col = mod_color))
  popViewport()
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(ggplot(bar_df, aes(x = Sample, y = Eigengene)) +
          geom_bar(stat = "identity", width = 1, fill = mod_color) +
          labs(y = "Eigengene", x = "Sample") + theme_bw(base_size = 14) +
          theme(panel.grid = element_blank(),
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
                axis.title = element_text(face = "bold"), plot.margin = margin(0, 5, 0, 5)),
        newpage = FALSE)
  popViewport()
  dev.off()
}

cat("All analyses completed.\n")