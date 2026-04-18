# 03_composition_analysis_and_figures.R
# Pre-Post LPS scRNAseq Data Analysis
# Step 3: Cell type annotation, compositional analysis, and figures
#
# Note: This script reads from the v5 Seurat objects produced by
# convert_seurat_v4_to_v5.R. Run that script first if starting from
# the raw objects. Marker gene tables referenced in annotations were
# produced by 02_marker_finding.R.

# Libraries ----
library(Seurat)
library(dplyr)
library(tibble)
library(openxlsx)
library(ggplot2)
library(speckle)
library(limma)

# Load data ----
integrated <- readRDS("integrated_all_samples_v5.RDS")
seu_mac    <- readRDS("macrophage_recluster_v5.RDS")

# Cell type annotation ----

## Overall cell type labels ----
# Fine-grained cell type labels per cluster, assigned based on marker gene expression
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

# Broad cell type labels used for visualizations and compositional analyses
integrated$cluster_lab <- case_when(
  integrated$seurat_clusters %in% c(0, 1, 2, 3, 5, 7, 8, 9, 10) ~ "AM",
  integrated$seurat_clusters %in% c(4, 6, 17)                    ~ "Lymphs",
  integrated$seurat_clusters %in% c(11, 13)                      ~ "Proliferative",
  integrated$seurat_clusters %in% c(12, 15, 16, 18)              ~ "DC",
  TRUE                                                            ~ "Epi Mast"
)

## Macrophage cell type labels ----
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

# Color palettes ----
# Defined here for use in both compositional plots and figures
cell_type_color <- c(
  "AM"            = '#7FC97F',
  "Lymphs"        = '#FDC086',
  "Proliferative" = '#ce9662',
  "DC"            = '#BEAED4',
  "Epi Mast"      = '#F0027F'
)

# Colors follow subtype grouping: greens = RAM, blues = Intermediate, reds = Recruited
cell_type_color_mac <- c(
  "RAM IGF1"          = "#006D2C",
  "RAM HES2"          = "#2CA25F",
  "RAM PLTP"          = "#66C2A4",
  "RAM lo"            = "#99D8C9",
  "RAM hi"            = "#CCECE6",
  "RAM old"           = "#EDF8FB",
  "Int Inflam"        = "#08519C",
  "Int Chol"          = "#3182BD",
  "Int IFN"           = "#6BAED6",
  "Int MT"            = "#9ECAE1",
  "Recr CCR2"         = "#A50F15",
  "Recr FOLR2"        = "#DE2D26",
  "Recr FOLR2 mature" = "#FB6A4A",
  "Recr CX3CR1"       = "#FC9272"
)

# Set factor levels for macrophage labels to control legend/axis order
seu_mac$cluster_lab <- factor(seu_mac$cluster_lab, levels = names(cell_type_color_mac))

# Cell type proportion tests ----

## Overall cell type proportions over time ----
# Exclude 96h timepoint — only 2 subjects (3156, 6386), insufficient for modeling
integrated_sub <- integrated[, integrated$time != 96]

# Ensure time is a factor with baseline (0) as reference
integrated_sub$time <- factor(integrated_sub$time, levels = c("0", "24", "48", "72", "120"))

# Get logit-transformed proportions
props <- getTransformedProps(
  clusters  = integrated_sub$cluster_lab,
  sample    = integrated_sub$library_id,
  transform = "logit"
)

# Build sample-level metadata aligned to props columns
# One row per sample, in the same order as columns of props$TransformedProps
meta <- integrated_sub@meta.data %>%
  select(library_id, subject, time, sex) %>%
  distinct() %>%
  arrange(match(library_id, colnames(props$TransformedProps)))

# Sanity check — must be TRUE before proceeding
stopifnot(all(meta$library_id == colnames(props$TransformedProps)))

# Design matrix
design <- model.matrix(~ time + sex, data = meta)

# Duplicate correlation to account for paired subject structure
dupcor <- duplicateCorrelation(
  props$TransformedProps,
  design,
  block = meta$subject
)
dupcor$consensus.correlation

# Fit model
fit <- lmFit(
  props$TransformedProps,
  design,
  block       = meta$subject,
  correlation = dupcor$consensus.correlation
)

# Step 1: F-test for any difference over time
# Test all time coefficients jointly per cell type
time_cols <- grep("^time", colnames(design), value = TRUE)

fit_f <- contrasts.fit(fit, coefficients = which(colnames(design) %in% time_cols))
fit_f <- eBayes(fit_f, robust = TRUE)

f_results <- topTableF(fit_f, number = Inf, sort.by = "F")

# Step 2: Pairwise contrasts vs baseline
# Only interpret if FDR-adjusted F-test is significant for that cell type
contrast_matrix <- makeContrasts(
  time24_vs_0  = time24,
  time48_vs_0  = time48,
  time72_vs_0  = time72,
  time120_vs_0 = time120,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2, robust = TRUE)

# Assemble results into long-format table
# FDR (adj.P.Val) is applied across cell types within each contrast by limma
f_res <- topTableF(fit_f, number = Inf, sort.by = "F") %>%
  rownames_to_column("Cell Type") %>%
  select(`Cell Type`, F, `F FDR` = adj.P.Val)

extract_contrast <- function(fit, coef, label) {
  topTable(fit, coef = coef, number = Inf, sort.by = "none") %>%
    rownames_to_column("Cell Type") %>%
    mutate(
      Comparison                = label,
      `logit(FC)`               = logFC,
      `t Statistic`             = t,
      `Unadjusted P`            = P.Value,
      `FDR (across cell types)` = adj.P.Val
    ) %>%
    select(`Cell Type`, Comparison, `logit(FC)`, `t Statistic`,
           `Unadjusted P`, `FDR (across cell types)`)
}

results_table <- bind_rows(
  extract_contrast(fit2, "time24_vs_0",  "24h vs baseline"),
  extract_contrast(fit2, "time48_vs_0",  "48h vs baseline"),
  extract_contrast(fit2, "time72_vs_0",  "72h vs baseline"),
  extract_contrast(fit2, "time120_vs_0", "120h vs baseline")
) %>%
  left_join(f_res, by = "Cell Type") %>%
  select(`Cell Type`, F, `F FDR`, Comparison,
         `logit(FC)`, `t Statistic`, `Unadjusted P`, `FDR (across cell types)`) %>%
  arrange(`Cell Type`, Comparison)

write.xlsx(results_table, 'overall_compositional_analysis_results.xlsx', rowNames = FALSE)

## Macrophage subtype proportions over time ----
# Exclude 96h timepoint — only 2 subjects, insufficient for modeling
seu_mac_sub <- seu_mac[, seu_mac$time != 96]

# Ensure time is a factor with baseline (0) as reference
seu_mac_sub$time <- factor(seu_mac_sub$time, levels = c("0", "24", "48", "72", "120"))

# Get logit-transformed proportions
props <- getTransformedProps(
  clusters  = seu_mac_sub$cluster_lab_broad,
  sample    = seu_mac_sub$library_id,
  transform = "logit"
)

# Build sample-level metadata aligned to props columns
meta <- seu_mac_sub@meta.data %>%
  select(library_id, subject, time, sex) %>%
  distinct() %>%
  arrange(match(library_id, colnames(props$TransformedProps)))

# Sanity check — must be TRUE before proceeding
stopifnot(all(meta$library_id == colnames(props$TransformedProps)))

# Design matrix
design <- model.matrix(~ time + sex, data = meta)

# Duplicate correlation to account for paired subject structure
dupcor <- duplicateCorrelation(
  props$TransformedProps,
  design,
  block = meta$subject
)
dupcor$consensus.correlation

# Fit model
fit <- lmFit(
  props$TransformedProps,
  design,
  block       = meta$subject,
  correlation = dupcor$consensus.correlation
)

# Step 1: F-test for any difference over time
time_cols <- grep("^time", colnames(design), value = TRUE)

fit_f <- contrasts.fit(fit, coefficients = which(colnames(design) %in% time_cols))
fit_f <- eBayes(fit_f, robust = TRUE)

f_results <- topTableF(fit_f, number = Inf, sort.by = "F")

# Step 2: Pairwise contrasts vs baseline
contrast_matrix <- makeContrasts(
  time24_vs_0  = time24,
  time48_vs_0  = time48,
  time72_vs_0  = time72,
  time120_vs_0 = time120,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2, robust = TRUE)

# Assemble long-format results table
f_res <- topTableF(fit_f, number = Inf, sort.by = "F") %>%
  rownames_to_column("Cell Type") %>%
  select(`Cell Type`, F, `F FDR` = adj.P.Val)

results_table <- bind_rows(
  extract_contrast(fit2, "time24_vs_0",  "24h vs baseline"),
  extract_contrast(fit2, "time48_vs_0",  "48h vs baseline"),
  extract_contrast(fit2, "time72_vs_0",  "72h vs baseline"),
  extract_contrast(fit2, "time120_vs_0", "120h vs baseline")
) %>%
  left_join(f_res, by = "Cell Type") %>%
  select(`Cell Type`, F, `F FDR`, Comparison,
         `logit(FC)`, `t Statistic`, `Unadjusted P`, `FDR (across cell types)`) %>%
  arrange(`Cell Type`, Comparison)

write.xlsx(results_table, 'mac_compositional_analysis_results.xlsx', rowNames = FALSE)

# Figures ----

## UMAPs ----
p1 <- DimPlot(integrated, reduction = 'umap', group.by = "cluster_lab", raster = FALSE) +
  xlab('UMAP 1') +
  ylab('UMAP 2') +
  scale_colour_manual(values = cell_type_color, name = "Cell Lineage") +
  ggtitle('All Cells') +
  theme(text = element_text(family = "Helvetica"))

ggsave('Figure_overall_umap.pdf', p1, device = 'pdf', width = 6, height = 4, units = 'in')

p1 <- DimPlot(seu_mac, reduction = 'umap', group.by = "cluster_lab", raster = FALSE) +
  xlab('UMAP 1') +
  ylab('UMAP 2') +
  scale_colour_manual(values = cell_type_color_mac, name = "Cell Lineage") +
  ggtitle('Macrophages') +
  theme(text = element_text(family = "Helvetica"))

ggsave('Figure_mac_umap.pdf', p1, device = 'pdf', width = 6, height = 4, units = 'in')

## Dot plots ----

# Order clusters by broad cell type grouping for visualization
integrated$cluster_dotplot <- factor(
  integrated$seurat_clusters,
  levels = rev(c("0",  "1",  "2",  "3",  "5",  "7",  "9",  "8",  "10",  # AM
                 "4",  "6",  "17",                                          # Lymphs
                 "11", "13",                                                # Proliferative
                 "12", "15", "16", "18",                                    # DC
                 "14"))                                                      # Epi Mast
)

# Canonical marker genes per broad cell type
markers <- c("MRC1", "HLA-DRA", "CD68", "SIGLEC1", "C5AR1",  # AM
             "CD3E", "LCK",                                     # Lymphs
             "TOP2A", "MKI67",                                  # Proliferative
             "CD1C", "BTLA",                                    # DC
             "MUC16", "CPA3")                                   # Epi Mast

# Normalize RNA assay for visualization (requires explicit layer specification in Seurat v5)
DefaultAssay(integrated) <- "RNA"
integrated <- NormalizeData(integrated, assay = "RNA", layer = "counts", verbose = FALSE)

p1 <- DotPlot(integrated, features = markers, group.by = "cluster_dotplot") +
  RotatedAxis() +
  ylab("Cluster") +
  theme(
    text        = element_text(family = "Helvetica"),
    axis.text.x = element_text(face = "italic")
  )

ggsave('Figure_overall_dotplot_rna.pdf', p1, device = 'pdf',
       width = 10, height = 4, units = 'in')

# Reverse factor level order for visualization (bottom to top on y-axis)
seu_mac$cluster_dotplot <- factor(seu_mac$cluster_lab,
                                  levels = rev(levels(seu_mac$cluster_lab)))

# Marker genes per macrophage subtype
markers_mac <- c(
  "INHBA",  "IGF1",   "HES2",   "FABP4",   # RAM IGF1, HES2
  "FTL",    "CD74",   "GLDN",   "HBEGF",   # RAM PLTP, lo
  "MALAT1", "NEAT1",  "CCL4",   "CXCL10",  # RAM hi, old
  "SQLE",   "LDLR",   "RSAD2",  "IFIT1",   # Int Chol, IFN
  "MT1M",   "MT1E",   "VCAN",   "FCN1",    # Int MT, Inflam
  "LGMN",   "MERTK",  "PLA2G7", "HIF1A",   # IMAM
  "CCR2",   "CX3CR1"                        # Recr CCR2, CX3CR1
)

# Normalize RNA assay for visualization (requires explicit layer specification in Seurat v5)
DefaultAssay(seu_mac) <- "RNA"
seu_mac <- NormalizeData(seu_mac, assay = "RNA", layer = "counts", verbose = FALSE)

p1 <- DotPlot(seu_mac, features = markers_mac, group.by = "cluster_dotplot") +
  RotatedAxis() +
  ylab("Cluster") +
  theme(
    text        = element_text(family = "Helvetica"),
    axis.text.x = element_text(face = "italic")
  )

ggsave('Figure_mac_dotplot_rna.pdf', p1, device = 'pdf',
       width = 10, height = 4, units = 'in')
