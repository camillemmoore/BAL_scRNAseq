# 02_marker_finding.R
# Pre-Post LPS scRNAseq Data Analysis
# Step 2: Marker gene identification for overall and macrophage clusters
#
# Note: This analysis was performed using Seurat v4. The resulting marker
# tables are provided as supplementary files. Users should have run
# convert_seurat_v4_to_v5.R before running this script.

# Libraries ----
library(Seurat)
library(openxlsx)

# Load data ----
# Read in v5 converted objects
integrated <- readRDS("integrated_all_samples_v5.RDS")
seu_mac    <- readRDS("macrophage_recluster_v5.RDS")

# Overall cluster marker finding ----
DefaultAssay(integrated) <- "SCT"
nClusters <- length(unique(integrated$seurat_clusters))

## Preliminary SCT markers ----
# Naive marker finding ignoring sample structure
combined <- list()
for (i in 1:nClusters) {
  temp <- FindMarkers(integrated, ident.1 = (i - 1),
                      min.pct = 0.1, logfc.threshold = 0.25,
                      only.pos = TRUE, assay = 'SCT')
  combined[[i]] <- temp[temp$p_val_adj < 0.05, ]
}
names(combined) <- paste0('cl_', seq(0, nClusters - 1))

write.xlsx(combined, 'joint_sct_combined_markers_list.xlsx', rowNames = TRUE)

## Conserved markers ----
# Markers tested within each sample
# Excludes samples with fewer than 3 cells in a given cluster
markers.list <- list()
sample_freq <- table(integrated$seurat_clusters, integrated$library_id)
lib.ids <- as.character(colnames(sample_freq))

for (i in 1:nClusters) {
  low_count <- sample_freq[i, ] < 3
  
  obj <- if (sum(low_count) == 0) {
    integrated
  } else {
    integrated[, !(integrated$library_id %in% lib.ids[low_count])]
  }
  
  temp <- FindConservedMarkers(obj,
                               ident.1 = i - 1,
                               min.pct = 0, logfc.threshold = 0, only.pos = FALSE,
                               assay = 'SCT',
                               features = rownames(combined[[i]]),
                               grouping.var = 'library_id')
  
  markers.list[[i]] <- temp
  
  # Count how many samples show significant adjusted p-value per marker
  if (nrow(markers.list[[i]]) > 0) {
    markers.list[[i]]$n_samples_significant <- rowSums(
      markers.list[[i]][, grep("p_val_adj", colnames(markers.list[[i]]))] < 0.05
    )
  }
}

names(markers.list) <- paste0('cl_', seq(0, nClusters - 1))
write.xlsx(markers.list, 'joint_sct_sample_specific_markers_list.xlsx', rowNames = TRUE)

# Macrophage cluster marker finding ----
DefaultAssay(seu_mac) <- "SCT"
nClusters <- length(unique(seu_mac$seurat_clusters))

## Preliminary SCT markers ----
combined <- list()
for (i in 1:nClusters) {
  temp <- FindMarkers(seu_mac, ident.1 = (i - 1),
                      min.pct = 0.1, logfc.threshold = 0.25,
                      only.pos = TRUE, assay = 'SCT')
  combined[[i]] <- temp[temp$p_val_adj < 0.05, ]
}
names(combined) <- paste0('cl_', seq(0, nClusters - 1))
write.xlsx(combined, 'mac_recluster_sct_combined_markers.xlsx', rowNames = TRUE)

## Conserved markers ----
# Excludes samples with fewer than 3 cells in a given cluster
markers.list <- list()
sample_freq <- table(seu_mac$seurat_clusters, seu_mac$library_id)
library_ids <- as.character(colnames(sample_freq))

for (i in 1:nClusters) {
  low_count <- sample_freq[i, ] < 3
  
  obj <- if (sum(low_count) == 0) {
    seu_mac
  } else {
    seu_mac[, !(seu_mac$library_id %in% library_ids[low_count])]
  }
  
  temp <- FindConservedMarkers(obj,
                               ident.1 = i - 1,
                               min.pct = 0, logfc.threshold = 0, only.pos = FALSE,
                               assay = 'SCT',
                               features = rownames(combined[[i]]),
                               grouping.var = 'library_id')
  
  markers.list[[i]] <- temp
  
  if (nrow(markers.list[[i]]) > 0) {
    markers.list[[i]]$n_samples_significant <- rowSums(
      markers.list[[i]][, grep("p_val_adj", colnames(markers.list[[i]]))] < 0.05
    )
  }
}

names(markers.list) <- paste0('cl_', seq(0, nClusters - 1))
write.xlsx(markers.list, 'mac_recluster_sct_sample_specific_markers.xlsx', rowNames = TRUE)
