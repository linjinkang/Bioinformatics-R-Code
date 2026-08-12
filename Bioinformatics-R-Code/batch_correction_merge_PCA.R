library(data.table)
library(limma)
library(ggplot2)
library(tools)
library(RColorBrewer)

output_dir <- file.path(".", "output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

file_list <- list.files(".", pattern = "\\.csv$", full.names = TRUE)
if (length(file_list) == 0) stop("No CSV files found!")

data_list <- lapply(file_list, function(f) {
  dt <- fread(f)
  setnames(dt, 1, "gene_symbol")
  return(dt)
})
merged_data <- Reduce(function(x, y) merge(x, y, by = "gene_symbol"), data_list)

expr_mat <- as.matrix(merged_data[, -1, with = FALSE])
rownames(expr_mat) <- merged_data$gene_symbol
mode(expr_mat) <- "numeric"

batch_ids <- file_path_sans_ext(basename(file_list))
sample_n <- sapply(data_list, function(x) ncol(x) - 1)
batch_vec <- unlist(mapply(rep, batch_ids, sample_n))
if (length(batch_vec) != ncol(expr_mat)) stop("Batch vector length mismatch!")

corrected_expr <- removeBatchEffect(expr_mat, batch = batch_vec)
corrected_df <- data.frame(gene_symbol = rownames(expr_mat), corrected_expr, check.names = FALSE)
write.csv(corrected_df, file = file.path(output_dir, "merged.csv"), row.names = FALSE)

plot_pca <- function(data_matrix, batch_vector, plot_title, out_file) {
  pca_res <- prcomp(t(data_matrix), scale. = TRUE)
  var_explained <- pca_res$sdev^2 / sum(pca_res$sdev^2)
  pca_df <- data.frame(
    sample_id = rownames(pca_res$x),
    pc1 = pca_res$x[, 1],
    pc2 = pca_res$x[, 2],
    batch = factor(batch_vector, levels = unique(batch_vector))
  )
  n_groups <- length(unique(batch_vector))
  palette_colors <- if (n_groups <= 8) brewer.pal(n_groups, "Set1") else colorRampPalette(brewer.pal(9, "Set1"))(n_groups)
  
  x_lab <- paste0("PC1 (", round(var_explained[1] * 100, 1), "%)")
  y_lab <- paste0("PC2 (", round(var_explained[2] * 100, 1), "%)")
  
  p <- ggplot(pca_df, aes(x = pc1, y = pc2, color = batch)) +
    geom_point(size = 2, alpha = 0.9) +
    stat_ellipse(aes(fill = batch), geom = "polygon", alpha = 0.14, linetype = 2, show.legend = FALSE) +
    scale_color_manual(values = palette_colors) +
    scale_fill_manual(values = palette_colors) +
    labs(title = plot_title, x = x_lab, y = y_lab, color = "Batch") +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.major = element_line(colour = "gray90", linetype = 2),
      legend.title = element_text(face = "plain"),
      legend.background = element_rect(colour = "black", fill = NA, size = 0.25),
      plot.title = element_text(size = 14, face = "plain"),
      axis.title = element_text(face = "plain")
    )
  ggsave(out_file, p, width = 4.5, height = 3.2)
}

plot_pca(expr_mat, batch_vec, "PCA Before Batch Effect Removal", file.path(output_dir, "PCA_Before.pdf"))
plot_pca(corrected_expr, batch_vec, "PCA After Batch Effect Removal", file.path(output_dir, "PCA_After.pdf"))
