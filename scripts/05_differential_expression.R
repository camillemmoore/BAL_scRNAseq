# 05_differential_expression.R
# Pre-Post LPS scRNAseq Data Analysis
# Run differential expression across time (within clusters) or between clusters
#


# Libraries ----
library(Seurat)
library(SingleCellExperiment)
library(REBEL)
library(DESeq2)
library(dplyr)
library(pheatmap)
library(enrichR)
library(edgeR)

# Read/format data ----
seu_mac    <- readRDS("macrophage_recluster_v5.RDS")

## Broad cluster labels
seu_mac$cluster_lab_broad <- case_when(
  seu_mac$seurat_clusters %in% c(0, 2, 3, 5, 12, 13) ~ "RAM",
  seu_mac$seurat_clusters %in% c(7, 9, 10, 11)        ~ "Intermediate",
  seu_mac$seurat_clusters %in% c(1, 8)                 ~ "MoAM",
  seu_mac$seurat_clusters %in% c(4, 6)                 ~ "IMAM"
)

## Convert to singlecellexperiment
sce_mac<-as.SingleCellExperiment(seu_mac, assay = "RNA")



# DE between timepoints----
## Run within each cluster

## Format data----

## Add day labels, combine days 4/5
sce_mac$day<-paste0("day_", sce_mac$time/24)
sce_mac$day<-factor(ifelse(sce_mac$day %in% c("day_4", "day_5"), "day_4_5", sce_mac$day),
                    levels=c(paste0("day_", 0:3), "day_4_5"))

## Run analysis for each cluster----

## contrasts for testing DE at each day vs. baseline
contrasts=list(
  day_1 = c(0,1,0,0,0,0),
  day_2 = c(0,0,1,0,0,0),
  day_3 = c(0,0,0,1,0,0),
  "day_4_5" = c(0,0,0,0,1,0)
) 

## Initialize lists to hold results
timepoint_de_test_res<-list()
timepoint_de_hms<-list()

## Loop through clusters
for(cluster in levels(seu_mac$seurat_clusters)){
  
  ### Filter/transform data ----
  ## Filter down to cluster of interest
  sce_mac_fil<-sce_mac[,sce_mac$seurat_clusters==cluster]
  
  ## Pseudobulk by sample (library_id)
  pb_dat<-rebelAggregate(sce_mac_fil, sampleVariable = "library_id")
  
  # Filter out genes without at least 1 CPM in 3 samples
  cpm_counts=cpm(pb_dat@assays@data$counts)
  cpm_more1<-apply(cpm_counts, 1, function(x) sum(x>1))
  idx_keep=which(cpm_more1>=3)
  pb_dat_fil=pb_dat[idx_keep,]
  
  # VST Transformation
  dds <- DESeqDataSet(pb_dat_fil, design=~day+sex)
  dds <- DESeq(dds)
  vsd.fixed <- varianceStabilizingTransformation(dds, blind=F)
  assays(pb_dat_fil)[["normcounts"]] <- assay(vsd.fixed)
  
  ### Fit models/run tests ----
  fitObj<-rebelFit(pb_dat_fil, fixedEffects = ~day+sex, subjectVariable = "subject",
                   pseudoBulk = T,parallel = F)
  
  ## Run joint test for differences between any day vs. baseline
  timepoint_de_test_res[[paste0("cluster_", cluster)]]<-rebelTest(fitObj, contrast=do.call(rbind, contrasts))
  
  ### Make heatmaps ----
  
  ## Filter VST data to significant genes
  genes_for_hm<-rownames(timepoint_de_test_res[[paste0("cluster_", cluster)]])[timepoint_de_test_res[[paste0("cluster_", cluster)]]$p_val_adj<.05]
  vst_fil=assay(vsd.fixed)[genes_for_hm, ]
  
  ## Get averages for each day
  avg_norm<-lapply(unique(pb_dat$day), function(day){
    vst_sub<-vst_fil[,pb_dat$day==day, drop=FALSE]
    rowMeans(vst_sub)
  })%>%bind_cols()%>%data.frame()
  colnames(avg_norm)<-unique(pb_dat$day)
  rownames(avg_norm)<-rownames(vst_fil)
  avg_norm<-avg_norm[,c("day_0", "day_1", "day_2", "day_3", "day_4_5")]
  
  ## Make heatmap - save to results
  timepoint_de_hms[[paste0("cluster_", cluster)]]=pheatmap(as.matrix(avg_norm),
                                                           show_rownames = F, scale="row", show_colnames = T,cluster_cols = F,
                                                           clustering_method="ward.D",  main = paste0("Cluster ", cluster))
  
  
}






# Between broad cluster DE ----

## Filter/pseudobulk/transform data ----

## Not interested in intermediate
sce_mac_fil<-sce_mac[,sce_mac$cluster_lab_broad!="Intermediate"]

## Save sample/broad cluster combined variable
sce_mac_fil$sample_broad_cluster<-paste0(sce_mac_fil$library_id, "_", 
                                         sce_mac_fil$cluster_lab_broad)


## Pseudobulk by sample/broad cluster
pb_dat<-rebelAggregate(sce_mac_fil, sampleVariable = "sample_broad_cluster")

# Filter out lowly expressed genes
cpm_counts=edgeR::cpm(pb_dat@assays@data$counts)
cpm_more1<-apply(cpm_counts, 1, function(x) sum(x>1))
idx_keep=which(cpm_more1>=min(table(pb_dat$cluster_lab_broad)))
pb_dat_fil=pb_dat[idx_keep,]

## VST transformation
dds <- DESeqDataSet(pb_dat_fil, design=~cluster_lab_broad+sex)
dds <- DESeq(dds)
vsd.fixed <- varianceStabilizingTransformation(dds, blind=F)
assays(pb_dat_fil)[["normcounts"]] <- assay(vsd.fixed)

## Fit models/run tests----
## Run model with random effects for sample (library id) and subject
fitObj<-rebelFit(pb_dat_fil, fixedEffects = ~cluster_lab_broad+sex, subjectVariable = "subject",
                 sampleVariable = "library_id",
                 pseudoBulk = T,parallel = F)

## Test all pairwise comparisons
contrasts=list(
  MoAM_v_IAM =  c(0,1,0,0),
  RAM_v_IAM =  c(0,0,1,0),
  RAM_v_MoAM = c(0,-1,1,0)
)

broad_cluster_de_res<-lapply(contrasts, function(contrast){
  test<-rebelTest(fitObj, contrast = contrast)%>%
    arrange(p_val_adj)
})

## Enrichment Analysis----

### Define enrichment function----

## df - results data.frame
## logFC_thresh - logFC threshold for significance
## pval_thresh - p-value threshold for significance
## p_val_type - column name for p-value (raw or adjusted)
## dir - direction of genes to test, up (upregulated) or down (downregulated)
## libraries - enrichment gene set libraries to test
find_enrichment=function(df, logFC_thresh, pval_thresh, p_val_type, dir,
                         libraries){
  df$my_pvals=df[,p_val_type]
  if(dir=="up"){
    genes=row.names(df%>%filter(my_pvals<pval_thresh, Estimate>logFC_thresh))
  }else genes=row.names(df%>%filter(my_pvals<pval_thresh, Estimate< -logFC_thresh))
  Sys.sleep(1)
  enrichr_raw=enrichr(genes, enrichr.libraries)
  enrichr_df<-lapply(names(enrichr_raw), function(db_name){
    
    db_res=enrichr_raw[[db_name]]
    if(nrow(db_res)>0){
      db_res=db_res%>%mutate(database=db_name)%>%filter(Adjusted.P.value < 0.05)%>%
        dplyr::select(-Old.P.value, -Old.Adjusted.P.value)
      return(db_res)
    }else return(NULL)
  })%>%bind_rows()
  rm(enrichr_raw)
  return(enrichr_df)
}

### Run enrichment on all comparisons----

## Ontologies of interest
enrichr.libraries <- c("GO_Biological_Process_2023","GO_Molecular_Function_2023",
                       "GO_Cellular_Component_2023",
                       "Reactome_2022")

## Loop through contrasts - save enrichment results
broad_cluster_enrich<-lapply(names(broad_cluster_de_res), function(comp_name){
  enrich_list=lapply(c("up", "down"), function(dir){
    find_enrichment(df=broad_cluster_de_res[[comp_name]], logFC_thresh =.25, 
                    pval_thresh = .05, p_val_type = "p_val_adj", dir=dir,
                    libraries = enrichr.libraries)%>%
      arrange(Adjusted.P.value)
  })
  names(enrich_list)=paste0(comp_name, "_", c("upreg", "downreg"))
  enrich_list
})
broad_cluster_enrich=unlist(broad_cluster_enrich, recursive=F)







# Between fine cluster DE ----

## Filter/pseudobulk/transform data ----

## Not interested in intermediate
sce_mac_fil<-sce_mac[,sce_mac$seurat_clusters %in% c(1, 4, 6,8)]
sce_mac_fil$seurat_clusters<-droplevels(sce_mac_fil$seurat_clusters)

## Save sample/cluster labels
sce_mac_fil$sample_cluster<-paste0(sce_mac_fil$library_id, "_c", sce_mac_fil$seurat_clusters)

## Pseudobulk by sample/cluster
pb_dat<-rebelAggregate(sce_mac_fil, sampleVariable = "sample_cluster")

# Filter out lowly expressed genes
cpm_counts=edgeR::cpm(pb_dat@assays@data$counts)
cpm_more1<-apply(cpm_counts, 1, function(x) sum(x>1))
idx_keep=which(cpm_more1>=min(table(pb_dat$seurat_clusters)))
pb_dat_fil=pb_dat[idx_keep,]

## VST transformation
dds <- DESeqDataSet(pb_dat_fil, design=~seurat_clusters+sex)
dds <- DESeq(dds)
vsd.fixed <- varianceStabilizingTransformation(dds, blind=F)
assays(pb_dat_fil)[["normcounts"]] <- assay(vsd.fixed)

## Fit modes/run tests----
## Run model with random effects for sample (library id) and subject
fitObj<-rebelFit(pb_dat_fil, fixedEffects = ~seurat_clusters+sex, 
                 subjectVariable = "subject",
                 sampleVariable = "library_id",
                 pseudoBulk = T,parallel = F)

## Test all pairwise comparisons
contrasts=list(
  cluster4_v_cluster1 =  c(0,1,0,0,0),
  cluster6_v_cluster1 =  c(0,0,1,0,0),
  cluster8_v_cluster1 = c(0,0,0,1,0),
  cluster4_v_cluster6 =  c(0,1,-1,0,0),
  cluster4_v_cluster8 =  c(0,1,0,-1,0),
  cluster6_v_cluster8 =  c(0,0,1,-1,0)
)

fine_cluster_de_res<-lapply(contrasts, function(contrast){
  test<-rebelTest(fitObj, contrast = contrast)%>%
    arrange(p_val_adj)
})

## Enrichment Analysis----

### Run enrichment on all comparisons----

## Using function/libraries defined in previous section

## Loop through contrasts - save enrichment results
fine_cluster_enrich<-lapply(names(fine_cluster_de_res), function(comp_name){
  enrich_list=lapply(c("up", "down"), function(dir){
    find_enrichment(df=fine_cluster_de_res[[comp_name]], logFC_thresh =.25, 
                    pval_thresh = .05, p_val_type = "p_val_adj", dir=dir,
                    libraries = enrichr.libraries)%>%
      arrange(Adjusted.P.value)
  })
  names(enrich_list)=paste0(comp_name, "_", c("upreg", "downreg"))
  enrich_list
})
fine_cluster_enrich=unlist(fine_cluster_enrich, recursive=F)



