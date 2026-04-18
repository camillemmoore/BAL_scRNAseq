# 01_integration_and_clustering.R
# Pre-Post LPS scRNAseq Data Analysis
# Step 1: Data loading, QC, SCTransform, RPCA integration, and clustering
#
# Note: This analysis was performed using Seurat v4. The resulting objects
# were subsequently converted to Seurat v5 format for public release —
# see convert_seurat_v4_to_v5.R. Users wishing to reproduce downstream
# analyses should start with 03_analysis_and_figures.R using the provided
# v5 objects.

# Libraries ----
library(Seurat)
library(dplyr)
library(Matrix)
library(openxlsx)

# Load sample metadata ----
# Update path to reflect your local directory
samples <- read.xlsx('joint_subjectID_key.xlsx')
locs <- samples$path

# QC, filtering, and SCTransform per sample ----
sDat <- list()
mito.cut <- 0.25

for (ii in 1:length(locs)) {
  library_id <- samples$Library_ID[ii]
  
  temp <- Read10X(paste0(locs[ii], "/outs/filtered_feature_bc_matrix"))
  sDat[[ii]] <- CreateSeuratObject(counts = temp, min.cells = 0, min.features = 0,
                                   project = library_id)
  
  # Calculate ribosomal and mitochondrial gene percentages
  ribo.genes <- grep(pattern = "^RPS|^RPL", x = rownames(sDat[[ii]]),
                     ignore.case = TRUE, value = TRUE)
  mito.genes <- grep(pattern = "^MT-|^MRPS|^MRPL", x = rownames(sDat[[ii]]),
                     ignore.case = TRUE, value = TRUE)
  
  percent.ribo <- Matrix::colSums(GetAssayData(sDat[[ii]], slot = 'counts')[ribo.genes, ]) /
    Matrix::colSums(GetAssayData(sDat[[ii]], slot = 'counts'))
  percent.mito <- Matrix::colSums(GetAssayData(sDat[[ii]], slot = 'counts')[mito.genes, ]) /
    Matrix::colSums(GetAssayData(sDat[[ii]], slot = 'counts'))
  
  sDat[[ii]][['percent.mito']] <- percent.mito
  sDat[[ii]][['percent.ribo']] <- percent.ribo
  
  # Add sample-level metadata
  sDat[[ii]]$library_id <- library_id
  sDat[[ii]]$subject    <- samples$Subject.ID[ii]
  sDat[[ii]]$sex        <- samples$Sex[ii]
  sDat[[ii]]$time       <- samples$Timepoint[ii]
  
  # Filter: exclude mito-high cells and likely doublets (top 2% by UMI count)
  feature.genes <- setdiff(rownames(sDat[[ii]]), c(mito.genes, ribo.genes))
  sDat[[ii]] <- subset(sDat[[ii]][feature.genes, ],
                       subset = percent.mito < mito.cut &
                         nCount_RNA < quantile(sDat[[ii]]$nCount_RNA, 0.98))
  
  # Filter: exclude low-quality cells (bottom 2% by feature count)
  sDat[[ii]] <- subset(sDat[[ii]],
                       subset = nFeature_RNA > quantile(sDat[[ii]]$nFeature_RNA, 0.02))
  
  # Normalize using SCTransform, regressing out sequencing depth and mito percentage
  sDat[[ii]] <- SCTransform(sDat[[ii]], verbose = FALSE,
                            vars.to.regress = c("nCount_RNA", "percent.mito"))
}

names(sDat) <- samples$Library_ID

saveRDS(sDat, 'joint_data_sctransform.rds')

# Integration (RPCA) ----
# Select reference samples to anchor integration
ref_samples <- which(names(sDat) %in% c("6386A", "1114B", "1987B", "4893A", "95498A", "67723B"))

options(future.globals.maxSize = 100 * 1000 * 1024^2)

integration.features <- SelectIntegrationFeatures(object.list = sDat, nfeatures = 3000)

sDat <- lapply(sDat, function(x) RunPCA(x, features = integration.features, verbose = FALSE))
sDat <- PrepSCTIntegration(object.list = sDat, anchor.features = integration.features,
                           verbose = FALSE)

anchors <- FindIntegrationAnchors(sDat, dims = 1:30, normalization.method = "SCT",
                                  anchor.features = integration.features, k.anchor = 20,
                                  reduction = "rpca", reference = ref_samples)

integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

# Clustering ----
DefaultAssay(integrated) <- "integrated"
integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
ElbowPlot(integrated, ndims = 30)

nPCAs <- 25
integrated <- FindNeighbors(integrated, dims = 1:nPCAs)
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:nPCAs, n.neighbors = 50)
integrated <- FindClusters(integrated, reduction = "pca", resolution = 0.4, algorithm = 3)

# Save ----
saveRDS(integrated, "integrated_all_samples.RDS")

# Macrophage subclustering ----
# Subset alveolar macrophage clusters from main object for reclustering
mac_clusters <- c(0, 1, 2, 3, 5, 7, 8, 9, 10)
seu_mac <- integrated[, integrated$seurat_clusters %in% mac_clusters]

# Retain overall cluster label before reclustering
seu_mac$overall_cluster <- seu_mac$seurat_clusters

# Recluster within macrophage subset using same parameters as main clustering
nPCAs <- 25
DefaultAssay(seu_mac) <- "integrated"
seu_mac <- FindNeighbors(seu_mac, dims = 1:nPCAs)
seu_mac <- RunUMAP(seu_mac, reduction = "pca", dims = 1:nPCAs, n.neighbors = 50)
seu_mac <- FindClusters(seu_mac, reduction = "pca", resolution = 0.4, algorithm = 3)

saveRDS(seu_mac, 'macrophage_recluster.RDS')

