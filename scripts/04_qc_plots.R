# QC Violin Plots ----

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

# Helper function: build and plot QC violins ----
plot_qc_violins <- function(metadata, group_var) {
  metadata %>%
    dplyr::select(all_of(group_var), nCount_RNA, nFeature_RNA, 
                  percent.mito, percent.ribo) %>%
    mutate(
      log10_UMI   = log10(nCount_RNA),
      log10_Genes = log10(nFeature_RNA)
    ) %>%
    dplyr::select(all_of(group_var), log10_UMI, log10_Genes, 
                  percent.mito, percent.ribo) %>%
    pivot_longer(cols = -all_of(group_var),
                 names_to  = "metric",
                 values_to = "value") %>%
    mutate(metric = factor(metric,
                           levels = c("log10_UMI", "log10_Genes", 
                                      "percent.mito", "percent.ribo"),
                           labels = c("log10(Total UMI per Cell)",
                                      "log10(Number of Genes per Cell)",
                                      "Mitochondrial Transcripts (%)",
                                      "Ribosomal Transcripts (%)"))) %>%
    ggplot(aes(x = .data[[group_var]], y = value, fill = .data[[group_var]])) +
    geom_violin(scale = "width", trim = TRUE) +
    facet_wrap(~ metric, scales = "free_y", ncol = 2) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      strip.text      = element_text(face = "bold"),
      axis.title.x    = element_blank()
    ) +
    ylab(NULL)
}

# Overall object QC ----
integrated <- readRDS("integrated_all_samples_v5.RDS")

## Build sample labels ----
integrated@meta.data <- integrated@meta.data %>%
  mutate(
    subject_numeric = factor(subject,
                             levels = unique(subject),
                             labels = paste0("S", seq_along(unique(subject)))),
    sample_label = factor(
      paste0(subject_numeric, "-", time, "h"),
      levels = {
        distinct(., subject_numeric, time) %>%
          arrange(subject_numeric, as.numeric(time)) %>%
          mutate(label = paste0(subject_numeric, "-", time, "h")) %>%
          pull(label)
      }
    )
  )

## QC by sample ----
plot_qc_violins(integrated@meta.data, "sample_label")
ggsave("Figure_qc_by_sample.pdf", device = "pdf", width = 10, height = 8, units = "in")

## QC by cluster ----
plot_qc_violins(integrated@meta.data %>%
                  mutate(seurat_clusters = factor(seurat_clusters)),
                "seurat_clusters")
ggsave("Figure_qc_by_cluster.pdf", device = "pdf", width = 10, height = 8, units = "in")

## QC by broad cell type ----
plot_qc_violins(integrated@meta.data, "cluster_lab")
ggsave("Figure_qc_by_celltype.pdf", device = "pdf", width = 10, height = 8, units = "in")

# Macrophage object QC ----
seu_mac <- readRDS("macrophage_recluster_v5.RDS")

## Build sample labels ----
seu_mac@meta.data <- seu_mac@meta.data %>%
  mutate(
    subject_numeric = factor(subject,
                             levels = unique(subject),
                             labels = paste0("S", seq_along(unique(subject)))),
    sample_label = factor(
      paste0(subject_numeric, "-", time, "h"),
      levels = {
        distinct(., subject_numeric, time) %>%
          arrange(subject_numeric, as.numeric(time)) %>%
          mutate(label = paste0(subject_numeric, "-", time, "h")) %>%
          pull(label)
      }
    )
  )

## QC by sample ----
plot_qc_violins(seu_mac@meta.data, "sample_label")
ggsave("Figure_mac_qc_by_sample.pdf", device = "pdf", width = 10, height = 8, units = "in")

## QC by cluster ----
plot_qc_violins(seu_mac@meta.data %>%
                  mutate(seurat_clusters = factor(seurat_clusters)),
                "seurat_clusters")
ggsave("Figure_mac_qc_by_cluster.pdf", device = "pdf", width = 10, height = 8, units = "in")

## QC by fine-grained cell type ----
plot_qc_violins(seu_mac@meta.data, "cluster_lab")
ggsave("Figure_mac_qc_by_celltype.pdf", device = "pdf", width = 10, height = 8, units = "in")

## QC by broad cell type ----
plot_qc_violins(seu_mac@meta.data, "cluster_lab_broad")
ggsave("Figure_mac_qc_by_celltype_broad.pdf", device = "pdf", width = 10, height = 8, units = "in")

