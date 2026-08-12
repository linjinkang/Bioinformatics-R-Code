logfc_threshold <- 1
adjp_threshold <- 0.05

library(limma)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(grid)
library(scales)

raw_expr <- read.table("TCGA-KIRC.csv", header = TRUE, sep = ",", check.names = FALSE, stringsAsFactors = FALSE)
gene_names <- raw_expr[, 1]
expr_values <- as.matrix(raw_expr[, -1])
rownames(expr_values) <- gene_names
expr_values <- apply(expr_values, 2, as.numeric)

expr_norm <- normalizeBetweenArrays(limma::avereps(log2(expr_values + 1)))
mode(expr_norm) <- "numeric"

sample_ids <- colnames(expr_norm)
tissue_code <- sapply(strsplit(sample_ids, "-"), `[`, 4)
sample_type <- ifelse(grepl("^11", tissue_code), "normal", ifelse(grepl("^01", tissue_code), "cancer", NA))
keep <- !is.na(sample_type)
sample_ids <- sample_ids[keep]
sample_type <- sample_type[keep]

ctrl_idx <- sample_type == "normal"
treat_idx <- sample_type == "cancer"
expr_combined <- cbind(expr_norm[, ctrl_idx, drop = FALSE], expr_norm[, treat_idx, drop = FALSE])
n_ctrl <- sum(ctrl_idx)
n_treat <- sum(treat_idx)

group_factor <- factor(c(rep("normal", n_ctrl), rep("cancer", n_treat)), levels = c("normal", "cancer"))
design_matrix <- model.matrix(~ 0 + group_factor)
colnames(design_matrix) <- c("normal", "cancer")
fit_lm <- lmFit(expr_combined, design_matrix)
fit_contrast <- contrasts.fit(fit_lm, makeContrasts(cancer - normal, levels = design_matrix))
fit_contrast <- eBayes(fit_contrast)
all_degs <- topTable(fit_contrast, adjust.method = "fdr", number = Inf)

sig_degs <- all_degs[abs(all_degs$logFC) > logfc_threshold & all_degs$adj.P.Val < adjp_threshold, , drop = FALSE]
sig_degs$Regulation <- ifelse(sig_degs$logFC > 0, "Up", "Down")
sig_degs$FoldChange <- 2^(sig_degs$logFC)
write.csv(cbind(Gene = rownames(sig_degs), sig_degs), "DE_significant_genes.csv", row.names = FALSE)
write.table(rownames(sig_degs), "TCGA.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

volcano_df <- all_degs
volcano_df$Gene <- rownames(volcano_df)
x_lim <- max(abs(volcano_df$logFC), na.rm = TRUE) * 1.1
y_max <- max(-log10(volcano_df$adj.P.Val), na.rm = TRUE) * 1.05

up_sig <- volcano_df[volcano_df$logFC > logfc_threshold & volcano_df$adj.P.Val < adjp_threshold, ]
up_sig <- up_sig[order(abs(up_sig$logFC), decreasing = TRUE), ]
up_labels <- head(up_sig$Gene, 6)
down_sig <- volcano_df[volcano_df$logFC < -logfc_threshold & volcano_df$adj.P.Val < adjp_threshold, ]
down_sig <- down_sig[order(abs(down_sig$logFC), decreasing = TRUE), ]
down_labels <- head(down_sig$Gene, 6)
volcano_df$label <- ifelse(volcano_df$Gene %in% c(up_labels, down_labels), volcano_df$Gene, "")

n_up <- sum(volcano_df$logFC > logfc_threshold & volcano_df$adj.P.Val < adjp_threshold)
n_down <- sum(volcano_df$logFC < -logfc_threshold & volcano_df$adj.P.Val < adjp_threshold)

volcano_plot <- ggplot(volcano_df, aes(logFC, -log10(adj.P.Val))) +
  geom_point(aes(color = logFC, size = -log10(adj.P.Val)), alpha = 0.85) +
  geom_point(data = subset(volcano_df, label != ""), aes(size = -log10(adj.P.Val)),
             shape = 21, fill = NA, color = "black", stroke = 0.5) +
  scale_color_gradientn(colors = c("#053061", "#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B", "#67001F"),
                        values = rescale(c(-6, -4, -2, 0, 2, 4, 6)), limits = c(-x_lim, x_lim), name = "log2FC") +
  scale_size_continuous(range = c(0.3, 3.5), name = "-log10(p_val)", breaks = pretty(c(0, y_max), n = 5)[-1]) +
  geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), linetype = "dashed", color = "grey50", linewidth = 0.6) +
  geom_hline(yintercept = -log10(adjp_threshold), linetype = "dashed", color = "grey50", linewidth = 0.6) +
  geom_text_repel(data = subset(volcano_df, label != ""), aes(label = label), size = 6,
                  max.overlaps = Inf, box.padding = 0.5, point.padding = 0.3,
                  segment.color = "grey40", segment.size = 0.3, fontface = "italic",
                  color = "black", show.legend = FALSE) +
  annotate("segment", x = -x_lim * 0.55, xend = -x_lim * 0.9, y = y_max * 0.97, yend = y_max * 0.97,
           arrow = arrow(length = unit(0.25, "cm"), type = "closed"), color = "#2B83BA", linewidth = 1.2) +
  annotate("text", x = -x_lim * 0.72, y = y_max * 0.97, label = paste0("Down (", n_down, ")"),
           color = "#2B83BA", size = 6, vjust = -0.8) +
  annotate("segment", x = x_lim * 0.55, xend = x_lim * 0.9, y = y_max * 0.97, yend = y_max * 0.97,
           arrow = arrow(length = unit(0.25, "cm"), type = "closed"), color = "#D7191C", linewidth = 1.2) +
  annotate("text", x = x_lim * 0.72, y = y_max * 0.97, label = paste0("Up (", n_up, ")"),
           color = "#D7191C", size = 6, vjust = -0.8) +
  annotate("text", x = x_lim * 0.98, y = -log10(adjp_threshold), label = paste0("p = ", adjp_threshold),
           hjust = 1, vjust = -0.5, size = 6, color = "black") +
  scale_x_continuous(limits = c(-x_lim, x_lim), expand = c(0.02, 0)) +
  scale_y_continuous(limits = c(0, y_max), expand = c(0.02, 0)) +
  labs(title = "Volcano Plot", x = "avg_log2FC", y = "-log10(p_val)") +
  theme_bw(base_size = 20) +
  theme(plot.title = element_text(face = "plain", hjust = 0.5, size = 20),
        axis.title = element_text(face = "plain", size = 20),
        axis.text = element_text(color = "black", size = 20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1.2),
        legend.position = "right",
        legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.5),
        legend.key = element_rect(fill = "white"),
        legend.title = element_text(face = "plain", size = 16),
        legend.text = element_text(size = 16),
        legend.spacing.y = unit(0.3, "cm"),
        plot.margin = margin(15, 15, 10, 10)) +
  guides(color = guide_colorbar(title = "avg_log2FC", title.position = "top", title.hjust = 0.5,
                                barwidth = 1.2, barheight = 10, frame.colour = "black", frame.linewidth = 0.5,
                                ticks.colour = "black", ticks.linewidth = 0.5, order = 1),
         size = guide_legend(title = "-log10(p_val)", title.position = "top", title.hjust = 0.5,
                             override.aes = list(alpha = 1, color = "grey30"), order = 2))
ggsave("DE_volcano_gradient.pdf", volcano_plot, width = 7, height = 6, dpi = 300)

if (nrow(sig_degs) > 0) {
  sig_up <- sig_degs[sig_degs$logFC > 0, ]
  sig_down <- sig_degs[sig_degs$logFC < 0, ]
  sig_up <- sig_up[order(sig_up$logFC, decreasing = TRUE), ]
  sig_down <- sig_down[order(sig_down$logFC, decreasing = FALSE), ]
  top_n <- 20
  heat_genes <- c(head(rownames(sig_up), top_n), head(rownames(sig_down), top_n))
  
  group_heat <- factor(c(rep("Normal", n_ctrl), rep("Cancer", n_treat)), levels = c("Normal", "Cancer"))
  sample_order <- order(group_heat, decreasing = FALSE)
  heat_matrix <- expr_combined[heat_genes, sample_order, drop = FALSE]
  
  annotation_col <- data.frame(Group = group_heat[sample_order])
  rownames(annotation_col) <- colnames(heat_matrix)
  ann_colors <- list(Group = c(Normal = "#2B83BA", Cancer = "#D7191C"))
  heat_colors <- colorRampPalette(c("#053061", "#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B", "#67001F"))(100)
  
  heatmap_obj <- pheatmap(heat_matrix, scale = "row", color = heat_colors,
                          cluster_rows = TRUE, cluster_cols = FALSE,
                          show_rownames = TRUE, show_colnames = FALSE,
                          annotation_col = annotation_col, annotation_colors = ann_colors,
                          border_color = NA, fontsize_row = 10, fontsize_col = 10,
                          main = paste0("Top ", top_n, " Up- and Down-regulated DEGs"), silent = TRUE)
  
  grob_table <- heatmap_obj$gtable
  row_idx <- which(grob_table$layout$name == "row_names")
  if (length(row_idx) > 0) grob_table$grobs[[row_idx]]$gp$font <- 3
  pdf("DE_heatmap_top50.pdf", width = 6, height = 6)
  grid.draw(grob_table)
  dev.off()
}
