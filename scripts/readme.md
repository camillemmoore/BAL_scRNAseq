## Requirements
### R packages
- Seurat v5
- SeuratObject v5
- speckle
- limma
- dplyr
- tibble
- ggplot2
- openxlsx
- Matrix
- tidyr

Install packages:
```r
install.packages(c("dplyr", "tibble", "ggplot2", "openxlsx", "Matrix", "tidyr"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("limma")

install.packages("remotes")
remotes::install_github("phipsonlab/speckle")

install.packages("Seurat")
```

## Usage
### Starting from raw data
Scripts should be run in order:
1. `01_integration_and_clustering.R` — processes raw 10X data, performs QC, 
   SCTransform normalization, RPCA integration, and clustering. 
   **Note:** This analysis was originally performed in Seurat v4.
2. `convert_seurat_v4_to_v5.R` — converts objects to Seurat v5 format and 
   cleans metadata. Run this before scripts 02-04.
3. `02_marker_finding.R` — identifies cluster marker genes.
4. `03_analysis_and_figures.R` — cell type annotation, compositional analysis, 
   and figure generation.
5. `04_qc_plots.R` — QC violin plots.

### Starting from processed objects
Download the v5 Seurat objects from GEO and run scripts 03 and 04 directly.
Update file paths in scripts to reflect your local directory.
