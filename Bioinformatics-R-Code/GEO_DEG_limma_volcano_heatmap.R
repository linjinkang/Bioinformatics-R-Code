setwd("C:\\OneDrive\\文档\\TCGA\\D-RGs\\04GEO差异分析")

threshold_logFC <- 1
threshold_adjP  <- 0.05

library(limma)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(scales)

file_expr <- "merged.csv"
tmp_head <- readLines(file_expr, 1)
sep <- ifelse(grepl(",", tmp_head), ",", "\t")

expr_raw <- read.table(file_expr, header = TRUE, sep = sep,
                       check.names = FALSE, stringsAsFactors = FALSE)
rownames(expr_raw) <- expr_raw[, 1]
expr_mat <- expr_raw[, -1, drop = FALSE]

expr_mat <- as.matrix(expr_mat)
expr_mat <- apply(expr_mat, 2, as.numeric)
rownames(expr_mat) <- rownames(expr_raw)
colnames(expr_mat) <- colnames(expr_raw)[-1]

cat("Expression matrix:", dim(expr_mat)[1], "genes x", dim(expr_mat)[2], "samples\n")

sample_names <- colnames(expr_mat)
group_info <- ifelse(grepl("_con$", sample_names, ignore.case = TRUE), "Control",
                     ifelse(grepl("_tre$", sample_names, ignore.case = TRUE), "Disease", "Unknown"))
if (any(group_info == "Unknown")) stop("Unknown group in sample names!")
num_ctrl <- sum(group_info == "Control")
num_treat <- sum(group_info == "Disease")
cat("Groups: Control =", num_ctrl, ", Disease =", num_treat, "\n")

group_labels <- factor(group_info, levels = c("Control", "Disease"))
design_mat <- model.matrix(~ 0 + group_labels)
colnames(design_mat) <- c("Control", "Disease")

fit <- lmFit(expr_mat, design_mat)
contrast_mat <- makeContrasts(Disease - Control, levels = design_mat)
fit2 <- contrasts.fit(fit, contrast_mat)
fit2 <- eBayes(fit2)

all_diff_results <- topTable(fit2, adjust.method = "fdr", number = Inf)
cat("DEGs with FDR < 0.05:", sum(all_diff_results$adj.P.Val < 0.05), "\n")

write.csv(cbind(Gene = rownames(all_diff_results), all_diff_results),
          file = "results.csv", row.names = FALSE)
cat("Saved: results.csv\n")

significant_DEGs <- all_diff_results[
  abs(all_diff_results$logFC) > threshold_logFC &
    all_diff_results$adj.P.Val < threshold_adjP, , drop = FALSE]

output_DEGs <- cbind(Gene = rownames(significant_DEGs), significant_DEGs)
output_DEGs$SE <- ifelse(as.numeric(output_DEGs$t) != 0,
                         abs(as.numeric(output_DEGs$logFC) / as.numeric(output_DEGs$t)),
                         NA)
output_DEGs$Regulation <- ifelse(as.numeric(output_DEGs$logFC) > 0, "Up", "Down")

desired_order <- c("Gene", "logFC", "Regulation", "SE", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
existing_cols <- intersect(desired_order, colnames(output_DEGs))
output_DEGs <- output_DEGs[, c(existing_cols, setdiff(colnames(output_DEGs), existing_cols))]

write.csv(output_DEGs, file = "significant_genes.csv", row.names = FALSE)
cat("Saved: significant_genes.csv\n")

write.table(data.frame(gene = rownames(significant_DEGs)), file = "GEO.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
cat("Saved: GEO.txt\n")

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
    values = rescale(c(-6, -4, -2, 0, 2, 4, 6)),
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

ggsave("volcano.pdf", gradient_volcano, width = 7, height = 6, dpi = 300)
cat("Saved: volcano.pdf\n")

if(nrow(significant_DEGs) > 0) {
  sig_up <- significant_DEGs[significant_DEGs$logFC > 0, ]
  sig_down <- significant_DEGs[significant_DEGs$logFC < 0, ]
  sig_up <- sig_up[order(sig_up$logFC, decreasing = TRUE), ]
  sig_down <- sig_down[order(sig_down$logFC, decreasing = FALSE), ]
  
  top_n <- 20
  up_genes   <- head(rownames(sig_up), top_n)
  down_genes <- head(rownames(sig_down), top_n)
  heatmap_genes <- c(up_genes, down_genes)
  
  heatmap_mat <- expr_mat[heatmap_genes, , drop = FALSE]
  sample_order <- order(group_labels, decreasing = FALSE)
  heatmap_mat <- heatmap_mat[, sample_order, drop = FALSE]
  group_ordered <- group_labels[sample_order]
  
  annotation_col <- data.frame(Group = group_ordered)
  rownames(annotation_col) <- colnames(heatmap_mat)
  ann_colors <- list(Group = c(Control = "#2B83BA", Disease = "#D7191C"))
  
  heat_colors <- colorRampPalette(
    c("#053061", "#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B", "#67001F")
  )(100)
  
  italic_labels <- parse(text = paste0("italic('", heatmap_genes, "')"))
  
  pheatmap(heatmap_mat,
           scale = "row",
           color = heat_colors,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           show_rownames = TRUE,
           show_colnames = FALSE,
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           labels_row = italic_labels,
           border_color = NA,
           fontsize_row = 10,
           fontsize_col = 10,
           main = paste0("Top ", top_n, " Up- and Down-regulated DEGs"),
           filename = "heatmap.pdf",
           width = 6,
           height = 6)
  cat("Saved: heatmap.pdf\n")
} else {
  cat("No significant DEGs found, skipping heatmap.\n")
}