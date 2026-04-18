# Clean the mould data for GEO and sharing
library(Seurat)
library(SeuratObject)
library(dplyr) 

# Convert original overall object to v5 ----
## Load and convert to Seurat v5 ----
# Load integrated object 
integrated <- readRDS("integrated_all_samples.RDS")
integrated <- UpdateSeuratObject(integrated)
integrated[["RNA"]]        <- as(integrated[["RNA"]], "Assay5")
integrated[["SCT"]]        <- as(integrated[["SCT"]], "Assay5")
integrated[["integrated"]] <- as(integrated[["integrated"]], "Assay5")


## Add cell type annotations ----
# ── Cell type map: cluster to fine-grained label ──────────────────────────────
celltype_map <- c(
  "0"  = "RAM",
  "1"  = "RecAM",
  "2"  = "RAM",
  "3"  = "RAM",
  "4"  = "T cell",
  "5"  = "RecAM",
  "6"  = "T cell",
  "7"  = "RAM",
  "8"  = "RecAM",
  "9"  = "RAM",
  "10" = "RAM",
  "11" = "Cycling",
  "12" = "DC",
  "13" = "Cycling",
  "14" = "Epi and Mast",
  "15" = "DC",
  "16" = "pDC",
  "17" = "B cell",
  "18" = "Migratory DC"
)

integrated$overall_celltype <- celltype_map[as.character(integrated$seurat_clusters)]

integrated$cluster_lab <- case_when(
  integrated$seurat_clusters %in% c(0, 1, 2, 3, 5, 7, 8, 9, 10) ~ "AM",
  integrated$seurat_clusters %in% c(4, 6, 17)                    ~ "Lymphs",
  integrated$seurat_clusters %in% c(11, 13)                      ~ "Proliferative",
  integrated$seurat_clusters %in% c(12, 15, 16, 18)              ~ "DC",
  TRUE                                                            ~ "Epi Mast"
)

## Normalize RNA for visualization ----
integrated <- NormalizeData(integrated, assay = "RNA", layer = "counts", verbose = FALSE)

## Save cleaned v5 object ----
saveRDS(integrated, "integrated_all_samples_v5.RDS")



# Convert macrophage subcluster object to v5 ----
## Load and convert to Seurat v5 ----
seu_mac <- readRDS("macrophage_recluster.RDS")
seu_mac <- UpdateSeuratObject(seu_mac)
seu_mac[["RNA"]]        <- as(seu_mac[["RNA"]], "Assay5")
seu_mac[["SCT"]]        <- as(seu_mac[["SCT"]], "Assay5")
seu_mac[["integrated"]] <- as(seu_mac[["integrated"]], "Assay5")


## Cell type annotations ----
# Broad macrophage subtype labels
# RAM = resident alveolar macrophages
# Intermediate = intermediate state between RAM and MoAM
# MoAM = monocyte-derived alveolar macrophages (recruited)
# IMAM = interstitial macrophages
seu_mac$cluster_lab_broad <- case_when(
  seu_mac$seurat_clusters %in% c(0, 2, 3, 5, 12, 13) ~ "RAM",
  seu_mac$seurat_clusters %in% c(7, 9, 10, 11)        ~ "Intermediate",
  seu_mac$seurat_clusters %in% c(1, 8)                 ~ "MoAM",
  seu_mac$seurat_clusters %in% c(4, 6)                 ~ "IMAM"
)

# Fine-grained macrophage subtype labels
seu_mac$cluster_lab <- case_when(
  seu_mac$seurat_clusters == 0  ~ "RAM IGF1",
  seu_mac$seurat_clusters == 2  ~ "RAM HES2",
  seu_mac$seurat_clusters == 3  ~ "RAM PLTP",
  seu_mac$seurat_clusters == 5  ~ "RAM lo",
  seu_mac$seurat_clusters == 12 ~ "RAM hi",
  seu_mac$seurat_clusters == 13 ~ "RAM old",
  seu_mac$seurat_clusters == 7  ~ "Int Inflam",
  seu_mac$seurat_clusters == 9  ~ "Int Chol",
  seu_mac$seurat_clusters == 10 ~ "Int IFN",
  seu_mac$seurat_clusters == 11 ~ "Int MT",
  seu_mac$seurat_clusters == 1  ~ "Recr CCR2",
  seu_mac$seurat_clusters == 4  ~ "Recr FOLR2",
  seu_mac$seurat_clusters == 6  ~ "Recr FOLR2 mature",
  seu_mac$seurat_clusters == 8  ~ "Recr CX3CR1"
)

## Normalize RNA for visualization ----
DefaultAssay(seu_mac) <- "RNA"
seu_mac <- NormalizeData(seu_mac, assay = "RNA", layer = "counts", verbose = FALSE)

## Verify cell type assignments ----
table(seu_mac$seurat_clusters, seu_mac$cluster_lab_broad)
table(seu_mac$seurat_clusters, seu_mac$cluster_lab)

## Save cleaned v5 object ----
saveRDS(seu_mac, "macrophage_recluster_v5.RDS")


