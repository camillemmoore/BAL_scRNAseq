

# Human airspace macrophage signatures are conserved during sterile lung injury and repair


## Overview
This repository contains code for the analysis of single-cell RNA sequencing 
(scRNAseq) data from bronchoalveolar lavage (BAL) cells collected before and 
after endobronchial endotoxin challenge in healthy adults. The study examines 
transcriptional programming of airspace macrophage (AM) subsets during 
resolution of acute lung inflammation and tissue repair.

Fifteen subjects underwent bronchoscopic lavage before and at a pre-assigned 
time point after endobronchial exposure to bacterial endotoxin (24h, 48h, 72h, 
or 120h post-challenge). Single-cell RNA sequencing was performed on BAL cells 
to characterize longitudinal changes in AM transcriptional programming.

## Citation
Manuscript currently under revision. Citation will be updated upon publication.

Moore CM, McManus SA, Wynn EA, McClendon JD, Janssen WJ* and Mould KJ*. Human airspace macrophage signatures are conserved during sterile lung injury and repair.


## Data availability
Raw and processed data are available on GEO:
- Baseline samples: GSE151928
- Post-LPS samples: GSE300946

Processed Seurat objects (v5) are provided directly on GEO:
- `integrated_all_samples_v5.RDS` — all cell types, 15 subjects x 2 timepoints
- `macrophage_recluster_v5.RDS` — macrophage subcluster object

## Repository structure

    ├── scripts/
    │   ├── 01_integration_and_clustering.R   # Data loading, QC, RPCA integration, clustering
    │   ├── 02_marker_finding.R               # Marker gene identification
    │   ├── 03_analysis_and_figures.R         # Cell type annotation, compositional analysis, figures
    │   ├── 04_qc_plots.R                     # QC violin plots
    │   ├── 05_differential_expression.R      # Differential expression testing
    │   ├── 06_trajectory_analysis.RMD        # Pseudo-time Cell Rank analysis
    │   └── convert_seurat_v4_to_v5.R         # Object conversion and metadata cleaning
    ├── data/
    │   └── joint_subjectID_key.xlsx          # Sample metadata key
    └── README.md

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
   cleans metadata. Run this before scripts 02-06.
3. `02_marker_finding.R` — identifies cluster marker genes.
4. `03_analysis_and_figures.R` — cell type annotation, compositional analysis, 
   and figure generation.
5. `04_qc_plots.R` — QC violin plots.
6. `05_differential_expression.R` — DEG testing between macrophage subtypes and timepoints.
7. `06_trajectory_analysis.RMD` - trajectory analysis in macrophages. 

### Starting from processed objects
Download the v5 Seurat objects from GEO and run scripts 03 to 06 directly.
Update file paths in scripts to reflect your local directory.

## Notes on data processing
- Clustering and integration (script 01) were performed using Seurat v4. 
  Objects were subsequently converted to Seurat v5 for public release 
  (see `convert_seurat_v4_to_v5.R`).
- Two samples (subjects 28015 and 33402) had library IDs corrected after 
  initial processing. The corrected IDs are reflected in all provided objects 
  and in `joint_subjectID_key.xlsx`. See the key file for details.
- The 96h timepoint (n=2 subjects) was excluded from compositional analyses 
  due to insufficient sample size for modeling.

## Correspondence
Camille M. Moore, Phd  
mooreca@njhealth.org  
Center for Genes, Environment and Health  
National Jewish Health  

## License
Code: MIT License  
Data: CC BY 4.0



