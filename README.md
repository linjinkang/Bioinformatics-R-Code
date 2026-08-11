# bioinformatics-R-scripts

This repository contains the **R scripts** used for the bioinformatics analyses in our published SCI papers.  
All scripts are collected here to serve as the permanent code reference for the **Data Availability** section, ensuring full computational reproducibility.

## 📁 Repository Structure

All R scripts are placed directly in the root directory.  
Each script is named with a self‑explanatory identifier (e.g., `TCGA_DEG_limma_volcano_heatmap.R`).  
Scripts are self‑contained; the header comment inside each file clearly states the associated paper, its DOI, and the specific analysis purpose.

## 🧬 Scripts

- **`TCGA_DEG_limma_volcano_heatmap.R`**  
  Differential expression analysis using limma, Nature‑style gradient volcano plot, and a heatmap of the top up‑ and down‑regulated genes (gene names in italics). Outputs a CSV of significant DEGs and a plain text gene list.

> More scripts will be added as papers are published. Each script’s file header contains the full citation information.

## 🚀 How to Use

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/bioinformatics-R-scripts.git
