setwd("C:")

threshold_logFC <- 1
threshold_adjP  <- 0.05

library(limma)
library(ggplot2)
library(ggrepel)
library(scales)

cat("[Step 1] Reading and preprocessing...\n")
file_expr <- "TCGA-KIRC.csv"
raw_expr <- read.table(file_expr, header = TRUE, sep = ",", check.names = FALSE, stringsAsFactors = FALSE)

gene_names <- raw_expr[, 1]
expr_values <- as.matrix(raw_expr[, -1])
rownames(expr_values) <- gene_names
expr_values <- apply(expr_values, 2, as.numeric)
rownames(expr_values) <- gene_names
expr_values <- log2(expr_values + 1)
expr_values <- limma::avereps(expr_values)
expr_values <- normalizeBetweenArrays(expr_values)
mode(expr_values) <- "numeric"
numeric_expr <- expr_values

sample_names <- colnames(numeric_expr)
sample_type_code <- sapply(strsplit(sample_names, "-"), function(x) x[4])
sample_type <- ifelse(grepl("^11", sample_type_code), "normal",
                      ifelse(grepl("^01", sample_type_code), "cancer", NA))

valid_samples <- !is.na(sample_type)
sample_names_valid <- sample_names[valid_samples]
sample_type_valid <- sample_type[valid_samples]
ctrl_samples <- sample_names_valid[sample_type_valid == "normal"]
treat_samples <- sample_names_valid[sample_type_valid == "cancer"]

data_ctrl <- numeric_expr[, ctrl_samples, drop = FALSE]
data_treat <- numeric_expr[, treat_samples, drop = FALSE]
combined_expr <- cbind(data_ctrl, data_treat)

num_ctrl <- ncol(data_ctrl)
num_treat <- ncol(data_treat)

group_labels <- c(rep("normal", num_ctrl), rep("cancer", num_treat))
design_mat <- model.matrix(~0 + factor(group_labels))
colnames(design_mat) <- levels(factor(group_labels))

fit_initial <- lmFit(combined_expr, design_mat)
contrast_mat <- makeContrasts(Comparison = cancer - normal, levels = design_mat)
fit_contrasted <- contrasts.fit(fit_initial, contrast_mat)
fit_contrasted <- eBayes(fit_contrasted)
all_diff_results <- topTable(fit_contrasted, adjust.method = "fdr", number = Inf)

cat("[Step 4] Filtering and exporting significant DEGs...\n")
significant_DEGs <- all_diff_results[abs(all_diff_results$logFC) > threshold_logFC &
                                       all_diff_results$adj.P.Val < threshold_adjP, ]

significant_DEGs$Regulation <- ifelse(significant_DEGs$logFC > 0, "Up", "Down")
significant_DEGs$FoldChange <- 2^(significant_DEGs$logFC)

output_DEGs <- cbind(Gene = rownames(significant_DEGs), significant_DEGs)
write.csv(output_DEGs, file = "DE_significant_genes.csv", row.names = FALSE)

write.table(rownames(significant_DEGs), file = "TCGA.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

cat(sprintf("Significant DEGs: %d (Up %d, Down %d)\n",
            nrow(significant_DEGs),
            sum(significant_DEGs$Regulation == "Up"),
            sum(significant_DEGs$Regulation == "Down")))

cat("[Step 5] Drawing volcano plot...\n")
volcano_data <- all_diff_results
volcano_data$Gene <- rownames(volcano_data)

x_limit <- max(abs(volcano_data$logFC), na.rm = TRUE) * 1.1
y_max <- max(-log10(volcano_data$adj.P.Val), na.rm = TRUE) * 1.05

up_sig <- volcano_data[volcano_data$logFC > threshold_logFC & volcano_data$adj.P.Val < threshold_adjP, ]
up_sig <- up_sig[order(abs(up_sig$logFC), decreasing = TRUE), ]
up_label <- head(up_sig$Gene, 6)

down_sig <- volcano_data[volcano_data$logFC < -threshold_logFC & volcano_data$adj.P.Val < threshold_adjP, ]
down_sig <- down_sig[order(abs(down_sig$logFC), decreasing = TRUE), ]
down_label <- head(down_sig$Gene, 6)

label_genes <- c(up_label, down_label)
volcano_data$label <- ifelse(volcano_data$Gene %in% label_genes, volcano_data$Gene, "")

up_count <- sum(volcano_data$logFC > threshold_logFC & volcano_data$adj.P.Val < threshold_adjP)
down_count <- sum(volcano_data$logFC < -threshold_logFC & volcano_data$adj.P.Val < threshold_adjP)

gradient_volcano <- ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = logFC, size = -log10(adj.P.Val)), alpha = 0.85) +
  geom_point(data = subset(volcano_data, label != ""),
             aes(size = -log10(adj.P.Val)),
             shape = 21, fill = NA, color = "black", stroke = 0.5) +
  scale_color_gradientn(
    colors = c("#053061", "#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B", "#67001F"),
    values = scales::rescale(c(-6, -4, -2, 0, 2, 4, 6)),
    limits = c(-x_limit, x_limit),
    name = "log2FC"
  ) +
  scale_size_continuous(
    range = c(0.3, 3.5),
    name = "-log10(p_val)",
    breaks = pretty(c(0, y_max), n = 5)[-1]
  ) +
  geom_vline(xintercept = c(-threshold_logFC, threshold_logFC),
             linetype = "dashed", color = "grey50", linewidth = 0.6) +
  geom_hline(yintercept = -log10(threshold_adjP),
             linetype = "dashed", color = "grey50", linewidth = 0.6) +
  geom_text_repel(
    data = subset(volcano_data, label != ""),
    aes(label = label),
    size = 6,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey40",
    segment.size = 0.3,
    fontface = "italic",
    color = "black",
    show.legend = FALSE
  ) +
  annotate("segment", x = -x_limit * 0.55, xend = -x_limit * 0.9,
           y = y_max * 0.97, yend = y_max * 0.97,
           arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
           color = "#2B83BA", linewidth = 1.2) +
  annotate("text", x = -x_limit * 0.72, y = y_max * 0.97,
           label = paste0("Down (", down_count, ")"), color = "#2B83BA",
           size = 6, vjust = -0.8) +
  annotate("segment", x = x_limit * 0.55, xend = x_limit * 0.9,
           y = y_max * 0.97, yend = y_max * 0.97,
           arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
           color = "#D7191C", linewidth = 1.2) +
  annotate("text", x = x_limit * 0.72, y = y_max * 0.97,
           label = paste0("Up (", up_count, ")"), color = "#D7191C",
           size = 6, vjust = -0.8) +
  annotate("text", x = x_limit * 0.98, y = -log10(threshold_adjP),
           label = paste0("p = ", threshold_adjP),
           hjust = 1, vjust = -0.5, size = 6, color = "black") +
  scale_x_continuous(limits = c(-x_limit, x_limit), expand = c(0.02, 0)) +
  scale_y_continuous(limits = c(0, y_max), expand = c(0.02, 0)) +
  labs(title = "Volcano Plot", x = "avg_log2FC", y = "-log10(p_val)") +
  theme_bw(base_size = 20) +
  theme(
    plot.title = element_text(face = "plain", hjust = 0.5, size = 20),
    axis.title = element_text(face = "plain", size = 20),
    axis.text = element_text(color = "black", size = 20),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 1.2),
    legend.position = "right",
    legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.5),
    legend.key = element_rect(fill = "white"),
    legend.title = element_text(face = "plain", size = 16),
    legend.text = element_text(size = 16),
    legend.spacing.y = unit(0.3, "cm"),
    plot.margin = margin(15, 15, 10, 10)
  ) +
  guides(
    color = guide_colorbar(
      title = "avg_log2FC",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = 1.2,
      barheight = 10,
      frame.colour = "black",
      frame.linewidth = 0.5,
      ticks.colour = "black",
      ticks.linewidth = 0.5,
      order = 1
    ),
    size = guide_legend(
      title = "-log10(p_val)",
      title.position = "top",
      title.hjust = 0.5,
      override.aes = list(alpha = 1, color = "grey30"),
      order = 2
    )
  )

ggsave("DE_volcano_gradient.pdf", gradient_volcano, width = 7, height = 6, dpi = 300)
cat("Gradient volcano plot saved: DE_volcano_gradient.pdf\n")
cat("All analyses completed!\n")

cat("[Step 6] Drawing heatmap for top 50 DEGs...\n")
library(pheatmap)
library(grid)

if(nrow(significant_DEGs) > 0) {
  sig_up   <- significant_DEGs[significant_DEGs$logFC > 0, ]
  sig_down <- significant_DEGs[significant_DEGs$logFC < 0, ]
  sig_up   <- sig_up[order(sig_up$logFC, decreasing = TRUE), ]
  sig_down <- sig_down[order(sig_down$logFC, decreasing = FALSE), ]
  
  top_n <- 20
  up_genes   <- head(rownames(sig_up), top_n)
  down_genes <- head(rownames(sig_down), top_n)
  
  heatmap_genes <- c(up_genes, down_genes)
  
  group_labels  <- factor(c(rep("Normal", num_ctrl), rep("Cancer", num_treat)),
                          levels = c("Normal", "Cancer"))
  sample_order  <- order(group_labels, decreasing = FALSE)
  heatmap_mat   <- combined_expr[heatmap_genes, sample_order, drop = FALSE]
  group_ordered <- group_labels[sample_order]
  
  annotation_col <- data.frame(Group = group_ordered)
  rownames(annotation_col) <- colnames(heatmap_mat)
  
  ann_colors <- list(Group = c(Normal = "#2B83BA", Cancer = "#D7191C"))
  
  heat_colors <- colorRampPalette(
    c("#053061", "#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B", "#67001F")
  )(100)
  
  heatmap_obj <- pheatmap(heatmap_mat,
                          scale = "row",
                          color = heat_colors,
                          cluster_rows = TRUE,
                          cluster_cols = FALSE,
                          show_rownames = TRUE,
                          show_colnames = F,
                          annotation_col = annotation_col,
                          annotation_colors = ann_colors,
                          border_color = NA,
                          fontsize_row = 10,
                          fontsize_col = 10,
                          main = paste0("Top ", top_n, " Up- and Down-regulated DEGs"),
                          silent = TRUE)
  
  g <- heatmap_obj$gtable
  row_idx <- which(g$layout$name == "row_names")
  if(length(row_idx) > 0) {
    g$grobs[[row_idx]]$gp$font <- 3
  }
  
  pdf("DE_heatmap_top50.pdf", width = 6, height = 6)
  grid.draw(g)
  dev.off()
  
  cat("Heatmap of top DEGs saved (italic gene names): DE_heatmap_top50.pdf\n")
} else {
  cat("No significant DEGs found, skipping heatmap.\n")
}