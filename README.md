# Network-Based Pathway Enrichment Analysis in Cervical Cancer

**Author:** Bhavana Chowdary Kothapalli  
**Institution:** University of Skövde, Sweden  
**Programme:** Master's in Bioinformatics (Second Cycle, 30 credits)  


---

## Overview

This repository contains the complete, reproducible bioinformatics workflow for the Master's thesis:

> *Network-Based Pathway Enrichment Analysis to Improve Molecular Insight in Cervical Cancer*

The study integrates conventional GO/KEGG enrichment with protein-protein interaction (PPI) network analysis applied to the GSE63514 cervical cancer microarray dataset (14 normal vs. 15 tumour samples; Affymetrix HG-U133 Plus 2.0).

---

## Repository Structure

```
cervical-cancer-network-analysis/
├── README.md
├── Dockerfile
├── scripts/
│   ├── 01_data_acquisition.R
│   ├── 02_quality_control.R
│   ├── 03_normalization.R
│   ├── 04_differential_expression.R
│   ├── 05_go_enrichment.R
│   ├── 06_kegg_enrichment.R
│   ├── 07_string_network.R
│   ├── 08_network_topology.R
│   ├── 09_hub_identification.R
│   ├── 10_module_detection.R
│   ├── 11_comparative_analysis.R
│   └── 12_visualization.R
├── outputs/
│   └── (generated figures and result tables saved here)
├── session_info.txt
└── master_analysis.Rmd
```

---

## How to Reproduce the Analysis

### Option 1: Using Docker (Recommended)

Docker ensures the exact same R version and package versions are used.

```bash
# 1. Pull the Docker image
docker pull bhavana/cervical-cancer-network:latest

# 2. Clone this repository
git clone https://github.com/YOUR_USERNAME/cervical-cancer-network-analysis.git
cd cervical-cancer-network-analysis

# 3. Run all scripts inside Docker
docker run --rm -v $(pwd):/workspace bhavana/cervical-cancer-network:latest \
  Rscript /workspace/scripts/01_data_acquisition.R
```

### Option 2: Running Locally in R

```r
# Install required packages first (see session_info.txt for exact versions)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "GEOquery", "affy", "limma", "preprocessCore",
  "hgu133plus2.db", "org.Hs.eg.db", "clusterProfiler",
  "GOSemSim", "STRINGdb", "igraph", "ggplot2",
  "pheatmap", "EnhancedVolcano", "ggrepel"
))
```

Then run scripts sequentially:

```bash
Rscript scripts/01_data_acquisition.R
Rscript scripts/02_quality_control.R
Rscript scripts/03_normalization.R
Rscript scripts/04_differential_expression.R
Rscript scripts/05_go_enrichment.R
Rscript scripts/06_kegg_enrichment.R
Rscript scripts/07_string_network.R
Rscript scripts/08_network_topology.R
Rscript scripts/09_hub_identification.R
Rscript scripts/10_module_detection.R
Rscript scripts/11_comparative_analysis.R
Rscript scripts/12_visualization.R
```

Or run the full analysis as a single R Markdown document:

```bash
Rscript -e "rmarkdown::render('master_analysis.Rmd')"
```

---

## Software Requirements

| Software | Version |
|----------|---------|
| R | 4.3.1 |
| Bioconductor | 3.17 |
| limma | 3.56.2 |
| clusterProfiler | 4.8.2 |
| STRINGdb | 2.12.0 |
| igraph | 1.5.1 |
| GOSemSim | 2.26.1 |
| preprocessCore | 1.62.1 |
| hgu133plus2.db | 3.13.0 |

Full session information is provided in `session_info.txt`.

---

## Data

The GSE63514 dataset is downloaded programmatically from NCBI GEO in script `01_data_acquisition.R`. No manual download is required. All data used is publicly available and de-identified.

- **GEO Accession:** GSE63514  
- **Platform:** Affymetrix Human Genome U133 Plus 2.0 Array (GPL570)  
- **Samples:** 14 normal cervical epithelium + 15 cervical squamous cell carcinoma  
- **Source:** den Boon et al. (2015), PNAS

---

## Outputs

All figures and result tables are saved to the `outputs/` directory:

- `outputs/figures/` — PCA plots, volcano plots, heatmaps, network graphs
- `outputs/tables/` — DEG tables, enrichment results, hub gene rankings

---

## Citation

If you use this workflow, please cite:

> Kothapalli, B. C. (2025). *Network-Based Pathway Enrichment Analysis to Improve Molecular Insight in Cervical Cancer*. Master's Thesis, University of Skövde.

---

## License

This project is licensed under the MIT License.
