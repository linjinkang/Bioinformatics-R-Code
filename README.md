# Bioinformatics-R-Code

This repository contains the **R scripts** used for the bioinformatics analyses in our published SCI papers.  
All scripts are collected here to serve as the permanent code reference for the **Data Availability** section, ensuring full computational reproducibility.

## 📁 Repository Structure

All R scripts are located in the **`Bioinformatics-R-Code/`** folder.  
Each script is named with a self‑explanatory identifier (e.g., `TCGA_DEG_limma_volcano_heatmap.R`).

## 🧬 Scripts

- **`TCGA_DEG_limma_volcano_heatmap.R`**  
  Differential expression analysis using limma, Nature‑style gradient volcano plot, and a heatmap of the top up‑ and down‑regulated genes (gene names in italics). Outputs a CSV of significant DEGs and a plain text gene list.

- **`GEO_DEG_limma_volcano_heatmap.R`**  
  Differential expression analysis for GEO datasets using limma, with automatic group recognition by sample name suffix (`_con` / `_tre`). Generates a gradient volcano plot and a heatmap of the top DEGs.

- **`batch_correction_merge_PCA.R`**  
  Merging of multiple expression matrices (CSV files) from different batches, followed by batch effect correction using `limma::removeBatchEffect` and generation of PCA plots before and after correction. The corrected merged expression matrix is saved as `merged.csv`.

- **`WGCNA_analysis.R`**  
  Weighted Gene Co‑expression Network Analysis (WGCNA) pipeline including soft‑threshold selection, module detection, module‑trait correlation, gene significance and module membership analysis. Outputs gene lists per module, module‑trait heatmaps, MM vs GS scatter plots, and module expression heatmaps with eigengene bar plots.

- **`Intersection_Venn_UpSet.R`**  
  Identification of gene intersections across multiple gene lists (`.txt` files). Generates a Venn diagram (via `ggvenn`), an UpSet plot (via `ComplexUpset`), global intersection gene list, and pairwise intersection tables.

> More scripts will be added as papers are published.

## 🚀 How to Use

1. Navigate to the `Bioinformatics-R-Code/` folder and open the R script(s) of interest in RStudio or any R environment.
2. Install the required packages (listed at the top of each script).
3. Obtain the necessary raw data from the sources mentioned in the corresponding paper’s “Data Availability” section.
4. Run the scripts following the order and instructions described in the paper’s Methods section.

## 📦 Required R Packages

Each script includes its own `library()` calls.  
Commonly used packages across the scripts include:  
`limma`, `ggplot2`, `ggrepel`, `pheatmap`, `scales`, `data.table`, `tools`, `RColorBrewer`

## 📝 Citation

If you use these scripts in your research, **please cite the corresponding original paper(s)**.  

## 📂 Data Availability

The scripts in this repository process publicly available or appropriately deposited datasets.  
**No raw data are stored here.**  
Please refer to the “Data Availability” section of the relevant paper for data accession numbers (GEO, TCGA, ArrayExpress, etc.).

## ⚖️ License

This repository is licensed under the **MIT License** – see the [LICENSE](LICENSE) file for details.  
You are free to use, modify, and distribute the code, provided you retain the original copyright notice.

## 📧 Contact

For questions about specific scripts, please contact [lin_jinkang@163.com] or open an issue in this repository.
