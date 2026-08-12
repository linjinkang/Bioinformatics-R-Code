logfc_threshold <- 1
adjp_threshold <- 0.05

library(limma)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(scales)

expr_file <- "expr_file.csv"
header_line <- readLines(expr_file, 1)
sep_char <- ifelse(grepl(",", header_line), ",", "\t")

raw_data <- read.table(expr_file, header = TRUE, sep = sep_char,
                       check.names = FALSE, stringsAsFactors = FALSE)
rownames(raw_data) <- raw_data[, 1]
expr_matrix <- as.matrix(raw_data[, -1, drop = FALSE])
expr_matrix <- apply(expr_matrix, 2, as.numeric)
rownames(expr_matrix) <- rownames(raw_data)
colnames(expr_matrix) <- colnames(raw_data)[-1]

sample_ids <- colnames(expr_matrix)
groups <- ifelse(grepl("_con$", sample_ids, ignore.case = TRUE), "Control",
                 ifelse(grepl("_tre$", sample_ids, ignore.case = TRUE), "Disease", "Unknown"))
if (any(groups == "Unknown")) stop("Unknown group in sample names!")

group_factor <- factor(groups, levels = c("Control", "Disease"))
design_matrix <- model.matrix(~ 0 + group_factor)
colnames(design_matrix) <- c("Control", "Disease")

fit_lm <- lmFit(expr_matrix, design_matrix)
contrast_matrix <- makeContrasts(Disease - Control, levels = design_matrix)
fit_contrast <- contrasts.fit(fit_lm, contrast_matrix)
fit_contrast <- eBayes(fit_contrast)

diff_results <- topTable(fit_contrast, adjust.method = "fdr", number = Inf)

write.csv(cbind(Gene = rownames(diff_results), diff_results),
          file = "results.csv", row.names = FALSE)

sig_genes <- diff_results[
  abs(diff_results$logFC) > logfc_threshold &
    diff_results$adj.P.Val < adjp_threshold, , drop = FALSE]

sig_output <- cbind(Gene = rownames(sig_genes), sig_genes)
sig_output$SE <- ifelse(as.numeric(sig_output$t) != 0,
                         abs(as.numeric(sig_output$logFC) / as.numeric(sig_output$t)),
                         NA)
sig_output$Regulation <- ifelse(as.numeric(sig_output$logFC) > 0, "Up", "Down")

col_order <- c("Gene", "logFC", "Regulation", "SE", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
existing_cols <- intersect(col_order, colnames(sig_output))
sig_output <- sig_output[, c(existing_cols, setdiff(colnames(sig_output), existing_cols))]

write.csv(sig_output, file = "significant_genes.csv", row.names = FALSE)

write.table(data.frame(gene = rownames(sig_genes)), file = "GEO.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

volcano_df <- diff_results
volcano_df$Gene <- rownames(volcano_df)

x_lim <- max(abs(volcano_df$logFC), na.rm = TRUE) * 1.1
y_max <- max(-log10(volcano_df$adj.P.Val), na.rm = TRUE) * 1.05

up_sig <- volcano_df[volcano_df$logFC > logfc_threshold & volcano_df$adj.P.Val < adjp_threshold, ]
up_sig <- up_sig[order(abs(up_sig$logFC), decreasing = TRUE), ]
up_labels <- head(up_sig$Gene, 6)

down_sig <- volcano_df[volcano_df$logFC < -logfc_threshold & volcano_df$adj.P.Val < adjp_threshold, ]
down_sig <- down_sig[order(abs(down_sig$logFC), decreasing = TRUE), ]
down_labels <- head(down_sig$Gene, 6)

label_genes_vec <- c(up_labels, down_labels)
volcano_df$label <- ifelse(volcano_df$Gene %in% label_genes_vec, volcano_df$Gene, "")

n_up <- sum(volcano_df$logFC > logfc_threshold & volcano_df$adj.P.Val < adjp_threshold)
n_down <- sum(volcano_df$logFC < -logfc_threshold & volcano_df$adj.P.Val < adjp_threshold)

volcano_plot <- ggplot(volcano_df, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = logFC, size = -log10(adj.P.Val)), alpha = 0.85) +
  geom_point(data = subset(volcano_df, label != ""),
             aes(size = -log10(adj.P.Val)),
             shape = 21, fill = NA, color = "black", stroke = 0.5) +
  scale_color_gradientn(
    colors = c("#053061", "#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B", "#67001F"),
    values = rescale(c(-6, -4, -2, 0, 2, 4, 6)),
    limits = c(-x_lim, x_lim),
    name = "log2FC"
  ) +
  scale_size_continuous(
    range = c(0.3, 3.5),
    name = "-log10(p_val)",
    breaks = pretty(c(0, y_max), n = 5)[-1]
  ) +
  geom_vline(xintercept = c(-logfc_threshold, logfc_threshold),
             linetype = "dashed", color = "grey50", linewidth = 0.6) +
  geom_hline(yintercept = -log10(adjp_threshold),
             linetype = "dashed", color = "grey50", linewidth = 0.6) +
  geom_text_repel(
    data = subset(volcano_df, label != ""),
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
  annotate("segment", x = -x_lim * 0.55, xend = -x_lim * 0.9,
           y = y_max * 0.97, yend = y_max * 0.97,
           arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
           color = "#2B83BA", linewidth = 1.2) +
  annotate("text", x = -x_lim * 0.72, y = y_max * 0.97,
           label = paste0("Down (", n_down, ")"), color = "#2B83BA",
           size = 6, vjust = -0.8) +
  annotate("segment", x = x_lim * 0.55, xend = x_lim * 0.9,
           y = y_max * 0.97, yend = y_max * 0.97,
           arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
           color = "#D7191C", linewidth = 1.2) +
  annotate("text", x = x_lim * 0.72, y = y_max * 0.97,
           label = paste0("Up (", n_up, ")"), color = "#D7191C",
           size = 6, vjust = -0.8) +
  annotate("text", x = x_lim * 0.98, y = -log10(adjp_threshold),
           label = paste0("p = ", adjp_threshold),
           hjust = 1, vjust = -0.5, size = 6, color = "black") +
  scale_x_continuous(limits = c(-x_lim, x_lim), expand = c(0.02, 0)) +
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

ggsave("volcano.pdf", volcano_plot, width = 7, height = 6, dpi = 300)

if (nrow(sig_genes) > 0) {
  sig_up <- sig_genes[sig_genes$logFC > 0, ]
  sig_down <- sig_genes[sig_genes$logFC < 0, ]
  sig_up <- sig_up[order(sig_up$logFC, decreasing = TRUE), ]
  sig_down <- sig_down[order(sig_down$logFC, decreasing = FALSE), ]
  
  top_n <- 20
  up_genes <- head(rownames(sig_up), top_n)
  down_genes <- head(rownames(sig_down), top_n)
  heat_genes <- c(up_genes, down_genes)
  
  heat_matrix <- expr_matrix[heat_genes, , drop = FALSE]
  sample_order <- order(group_factor, decreasing = FALSE)
  heat_matrix <- heat_matrix[, sample_order, drop = FALSE]
  group_ordered <- group_factor[sample_order]
  
  annotation_col <- data.frame(Group = group_ordered)
  rownames(annotation_col) <- colnames(heat_matrix)
  ann_colors <- list(Group = c(Control = "#2B83BA", Disease = "#D7191C"))
  
  heat_colors <- colorRampPalette(
    c("#053061", "#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B", "#67001F")
  )(100)
  
  italic_labels <- parse(text = paste0("italic('", heat_genes, "')"))
  
  pheatmap(heat_matrix,
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
}
