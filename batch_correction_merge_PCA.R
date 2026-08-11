library(data.table)
library(limma)
library(ggplot2)
library(tools)
library(RColorBrewer)

data_path  <- "C:\\OneDrive\\文档\\TCGA\\D-RGs\\03GEO数据合并"
output_dir <- file.path(data_path, "output")
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
cat("Output directory:", output_dir, "\n")

file_list <- list.files(data_path, pattern = "\\.csv$", full.names = TRUE)
if (length(file_list) == 0) stop("No CSV files found!")

cat("Merging", length(file_list), "files...\n")
data_list <- lapply(file_list, function(f) {
  dt <- fread(f)
  setnames(dt, 1, "geneSymbol")
  return(dt)
})
merged_data <- Reduce(function(x, y) merge(x, y, by = "geneSymbol"), data_list)

expr_matrix <- as.matrix(merged_data[, -1, with = FALSE])
rownames(expr_matrix) <- merged_data$geneSymbol
mode(expr_matrix) <- "numeric"

batch_names  <- file_path_sans_ext(basename(file_list))
sample_counts <- sapply(data_list, function(x) ncol(x) - 1)
batch_info    <- unlist(mapply(rep, batch_names, sample_counts))

if (length(batch_info) != ncol(expr_matrix)) stop("Batch vector length mismatch!")

cat("Performing batch correction...\n")
corrected_expr <- removeBatchEffect(expr_matrix, batch = batch_info)

corrected_df <- data.frame(geneSymbol = rownames(expr_matrix), corrected_expr, check.names = FALSE)
write.csv(corrected_df, file = file.path(output_dir, "merged.csv"), row.names = FALSE)
cat("Saved: merged.csv\n")

pca_advanced_plot <- function(mat, batch_info, title, file) {
  pca_res <- prcomp(t(mat), scale. = TRUE)
  pc_var  <- pca_res$sdev^2 / sum(pca_res$sdev^2)
  pca_df  <- data.frame(
    Sample = rownames(pca_res$x),
    PC1    = pca_res$x[,1],
    PC2    = pca_res$x[,2],
    Batch  = factor(batch_info, levels = unique(batch_info))
  )
  n_batch <- length(unique(batch_info))
  mycol   <- if (n_batch <= 8) brewer.pal(n_batch, "Set1") else colorRampPalette(brewer.pal(9,"Set1"))(n_batch)
  
  pc1_lab <- paste0("PC1 (", round(pc_var[1]*100, 1),"%)")
  pc2_lab <- paste0("PC2 (", round(pc_var[2]*100, 1),"%)")
  
  p <- ggplot(pca_df, aes(PC1, PC2, color = Batch)) +
    geom_point(size = 2, alpha = 0.9) +
    stat_ellipse(aes(fill = Batch), geom = "polygon", alpha = 0.14, linetype = 2, show.legend = FALSE) +
    scale_color_manual(values = mycol) +
    scale_fill_manual(values = mycol) +
    labs(title = title, x = pc1_lab, y = pc2_lab, color = "Batch") +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.major = element_line(colour = "gray90", linetype = 2),
      legend.title     = element_text(face = "plain"),
      legend.background = element_rect(colour = "black", fill = NA, size = 0.25),
      plot.title       = element_text(size = 14, face = "plain"),
      axis.title       = element_text(face = "plain")
    )
  ggsave(file, p, width = 4.5, height = 3.2)
  invisible(p)
}

pca_advanced_plot(expr_matrix, batch_info,
                  "PCA Before Batch Effect Removal",
                  file.path(output_dir, "PCA_Before.pdf"))

pca_advanced_plot(corrected_expr, batch_info,
                  "PCA After Batch Effect Removal",
                  file.path(output_dir, "PCA_After.pdf"))

cat("PCA plots saved.\nAll steps completed!\n")