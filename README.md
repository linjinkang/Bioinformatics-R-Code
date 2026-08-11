# bioinformatics-R-scripts

This repository contains the **R scripts** used for the bioinformatics analyses in our published SCI papers.  
All scripts are collected here to serve as the permanent code reference for the **Data Availability** section, ensuring full computational reproducibility.

## 📁 Repository Structure

To keep things simple, all R scripts are placed directly in the root directory.  
Each script is named with a short identifier linking it to the corresponding paper (e.g., `2024_BreastCancer_DEG.R`).  
If a script is used across multiple papers, it is clearly noted in its header comment.

## 📜 List of Scripts & Associated Papers

| Script file | Paper (short reference) | Key analysis |
|-------------|------------------------|--------------|
| `2024_BreastCancer_DEG.R` | Breast cancer prognosis model (2024) | Differential expression, survival analysis |
| `2023_scRNAseq_clustering.R` | Single-cell atlas of lung cancer (2023) | Seurat clustering, marker identification |
| `2022_TCGA_methylation.R` | Pan-cancer methylation study (2022) | Methylation 450K array preprocessing |
| *(add your own rows)* | | |

> **Note:** If a script belongs to multiple papers, it is listed under each relevant paper with a note in the file header.

## 🚀 How to Use

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/bioinformatics-R-scripts.git
