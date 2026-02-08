################################################################################
# Immune Cell Profiling in Anti-GAD65 - Data Analysis and Visualization
################################################################################
#
# Author: Sumanta Barman
# Description: Complete R analysis pipeline for immune cell profiling in 
#              anti-GAD65 encephalitis including single-cell RNA-seq processing,
#              cell type annotation, differential abundance, differential 
#              expression, and immune repertoire analysis.
#
# GitHub: [Your repository URL]
# License: [Your license]
#
################################################################################

################################################################################
# TABLE OF CONTENTS
################################################################################
# 1. Setup and Library Loading
# 2. Data Loading
# 3. Exploratory Visualization
# 4. Cell Type Annotation
#    4.1. SingleR Annotation
#    4.2. Azimuth Annotation
#    4.3. Marker Gene Identification
#    4.4. Manual Annotation
# 5. Differential Cell Abundance Analysis
#    5.1. CSF Compartment
#    5.2. PBMC Compartment
# 6. Differential Gene Expression
# 7. Pseudobulk Volcano Plot
# 8. Session Information
################################################################################

################################################################################
# 1. SETUP AND LIBRARY LOADING
################################################################################

# Set global options
options(
  stringsAsFactors = FALSE,
  warn = -1  # Suppress warnings
)

# Single-cell analysis packages
library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(DoubletFinder)
library(reticulate)
library(SeuratData)
library(SeuratWrappers)
library(patchwork)
library(harmony)
library(scuttle)
library(sctransform)
library(SingleR)
library(celldex)
library(pheatmap)
library(Azimuth)
library(Signac)
library(dplyr)
library(cowplot)
library(enrichR)
library(rafalib)
library(Matrix)
library(edgeR)

# Differential abundance packages
library(speckle)
library(limma)

# Immune repertoire analysis
library(immunarch)

# Visualization packages
library(EnhancedVolcano)
library(ggrepel)
library(RColorBrewer)
library(ggsci)
library(ggpubr)
library(rstatix)

# Set language
Sys.setenv(lang = "en_US")

cat("Libraries loaded successfully!\n\n")

################################################################################
# 2. DATA LOADING
################################################################################

# Load the integrated Seurat object with singlets
cat("Loading preprocessed Seurat object...\n")
seurat_object <- readRDS(file = "seurat_object_combined_singlets_integrated.rds")

# Display object summary
cat("\nSeurat object loaded:\n")
print(seurat_object)
cat("\n")

################################################################################
# 3. EXPLORATORY VISUALIZATION
################################################################################

cat("Generating UMAP visualizations...\n")

# Visualize different metadata variables on UMAP
DimPlot(
  seurat_object,
  reduction = "umap.harmony",
  group.by = c("sample", "Center", "Disease", "harmony_clusters"),
  combine = FALSE,
  label.size = 2,
  label = TRUE
)

################################################################################
# 4. CELL TYPE ANNOTATION
################################################################################

# ---- 4.1. SingleR Annotation ----
cat("Running SingleR annotation...\n")

# Load reference dataset
ref <- celldex::DatabaseImmuneCellExpressionData()

# Prepare Seurat object for SingleR
seurat_object_singler <- seurat_object
seurat_object_singler.counts <- GetAssayData(
  seurat_object_singler,
  slot = "counts",
  assay = "SCT"
)

# Run SingleR prediction
pred <- SingleR(
  test = seurat_object_singler.counts,
  ref = ref,
  labels = ref$label.fine
)

# Add SingleR labels to metadata
seurat_object_singler$singleR.labels <- pred$labels[
  match(rownames(seurat_object_singler@meta.data), rownames(pred))
]

# Save annotated object
saveRDS(seurat_object_singler, file = "immune_cell_profiling_anti_gad65_singler.rds")
cat("SingleR annotation complete. Object saved.\n\n")

## SingleR Visualization
cat("Visualizing SingleR results...\n")

# Basic UMAP with SingleR labels
p1 <- DimPlot(
  seurat_object_singler,
  reduction = "umap.harmony",
  group.by = "singleR.labels"
)

# UMAP with labels and repelling
p2 <- DimPlot(
  seurat_object_singler,
  reduction = "umap.harmony",
  group.by = "singleR.labels",
  label = TRUE,
  repel = TRUE,
  label.size = 2
)

# Split by disease condition
p3 <- DimPlot(
  seurat_object_singler,
  reduction = "umap.harmony",
  group.by = "singleR.labels",
  split.by = "Disease",
  label = TRUE,
  repel = TRUE,
  label.size = 3
)

print(p1 | p2)
print(p3)

# ---- 4.2. Azimuth Annotation ----
cat("\nRunning Azimuth annotation...\n")

# Prepare object for Azimuth
seurat_object_azimuth <- seurat_object
seurat_object_azimuth <- PrepSCTFindMarkers(
  seurat_object_azimuth,
  assay = "SCT",
  verbose = TRUE
)

# Run Azimuth annotation
seurat_object_azimuth <- RunAzimuth(
  seurat_object_azimuth,
  reference = "pbmcref"
)

# Save annotated object
saveRDS(seurat_object_azimuth, file = "immune_cell_profiling_anti_gad65_azimuth.rds")
cat("Azimuth annotation complete. Object saved.\n\n")

## Azimuth Visualization
cat("Visualizing Azimuth results...\n")

# UMAP without legend
p1 <- DimPlot(
  seurat_object_azimuth,
  reduction = "umap.harmony",
  group.by = "predicted.celltype.l2",
  label = TRUE,
  repel = TRUE,
  label.size = 2
) + NoLegend()

# UMAP with legend
p2 <- DimPlot(
  seurat_object_azimuth,
  reduction = "umap.harmony",
  group.by = "predicted.celltype.l2",
  label = TRUE,
  repel = TRUE,
  label.size = 3
)

# Split by disease
p3 <- DimPlot(
  seurat_object_azimuth,
  reduction = "umap.harmony",
  group.by = "predicted.celltype.l2",
  split.by = "Disease",
  label = TRUE,
  repel = TRUE,
  label.size = 3
) + NoLegend()

print(p1 | p2)
print(p3)

# ---- 4.3. Marker Gene Identification ----
cat("\nIdentifying cluster marker genes...\n")

# Prepare for marker finding
seurat_object <- PrepSCTFindMarkers(seurat_object, assay = "SCT", verbose = TRUE)
DefaultAssay(seurat_object)

# Set identity to harmony clusters
Idents(seurat_object) <- "harmony_clusters"

# Find markers using Wilcoxon test
markers_genes_RNA_wilcox <- FindAllMarkers(
  seurat_object,
  test.use = "wilcox",
  only.pos = TRUE,
  assay = "SCT"
)

# Save markers
write.csv(
  markers_genes_RNA_wilcox,
  "cluster_markers_wilcox.csv",
  row.names = FALSE
)
cat("Marker genes saved to cluster_markers_wilcox.csv\n\n")

# ---- 4.4. Manual Annotation ----
cat("Applying manual cell type annotations...\n")

# Create annotated object
seurat_object_annotated <- seurat_object

# Define cell type labels for each cluster
replacement_names <- c(
  "0" = "CD4+ Naive T",
  "1" = "CD4+ TCM",
  "2" = "CD8+ TEM, MAIT, gDT",
  "3" = "NK1",
  "4" = "CD8+ TEM",
  "5" = "Granulo",
  "6" = "cDC2",
  "7" = "Activated CD4+ TSCM",
  "8" = "CD8+ Naive T",
  "9" = "cDC1",
  "10" = "BC1",
  "12" = "BC2",
  "13" = "CD16+ Mono",
  "14" = "NK2",
  "16" = "Mega",
  "17" = "pDC",
  "18" = "CD14+ Mono",
  "19" = "Memory BC",
  "20" = "Plasma",
  "21" = "Microglia-like",
  "24" = "BC3"
)

# Add annotations to metadata
metadata <- seurat_object_annotated@meta.data
metadata$annotate <- replacement_names[as.character(metadata$harmony_clusters)]
seurat_object_annotated@meta.data <- metadata

# Save annotated object
saveRDS(seurat_object_annotated, file = "immune_cell_profiling_anti_gad65_annotated.rds")
cat("Manual annotation complete. Object saved.\n\n")

## Visualization of Annotated Clusters
cat("Visualizing annotated cell types...\n")

# UMAP with cell type labels
p1 <- DimPlot(
  seurat_object_annotated,
  reduction = "umap.harmony",
  group.by = "annotate",
  label = TRUE
)

# Split by disease
p2 <- DimPlot(
  seurat_object_annotated,
  reduction = "umap.harmony",
  group.by = "annotate",
  split.by = "Disease"
)

# Split by sample type
p3 <- DimPlot(
  seurat_object_annotated,
  reduction = "umap.harmony",
  group.by = "annotate",
  split.by = "Sampletype"
)

print(p1)
print(p2)
print(p3)

################################################################################
# 5. DIFFERENTIAL CELL ABUNDANCE ANALYSIS
################################################################################

cat("\n=== DIFFERENTIAL CELL ABUNDANCE ANALYSIS ===\n\n")

# Display cell counts
cat("Cell type distribution:\n")
print(table(seurat_object_annotated$annotate))
cat("\nSample distribution:\n")
print(table(seurat_object_annotated$sample))

# ---- Proportions by Sample Type (CSF vs PBMC) ----
cat("\nAnalyzing proportions by sample type...\n")

props_csf_pbmc <- getTransformedProps(
  clusters = seurat_object_annotated$annotate,
  sample = seurat_object_annotated$Sampletype
)

# Visualize proportions
p1 <- plotCellTypeProps(
  clusters = seurat_object_annotated$annotate,
  sample = seurat_object_annotated$Sampletype
) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  ggtitle("Cell Type Proportions by Sample Type")

print(p1)

# ---- Proportions by Disease Status ----
cat("\nAnalyzing proportions by disease status...\n")

props_disease <- getTransformedProps(
  clusters = seurat_object_annotated$annotate,
  sample = seurat_object_annotated$Disease
)

# Visualize
p1 <- plotCellTypeProps(
  clusters = seurat_object_annotated$annotate,
  sample = seurat_object_annotated$Disease
) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  ggtitle("Cell Type Proportions by Disease Status")

print(p1)

# ---- Subset CSF and PBMC Compartments ----
cat("\nSubsetting CSF and PBMC compartments...\n")

CSF <- subset(seurat_object_annotated, Sampletype == "CSF")
cat(sprintf("CSF samples: %d cells\n", ncol(CSF)))

PBMC <- subset(seurat_object_annotated, Sampletype == "PBMC")
cat(sprintf("PBMC samples: %d cells\n\n", ncol(PBMC)))

# ---- 5.1. CSF Compartment Analysis ----
cat("=== CSF COMPARTMENT ANALYSIS ===\n")

## Cell Type Proportions in CSF
props_csf <- getTransformedProps(
  clusters = CSF$annotate,
  sample = CSF$samplename
)

# Save proportions
write.csv(
  props_csf,
  "cluster_abundances_csf_annotated.csv",
  quote = FALSE,
  row.names = FALSE
)
cat("CSF proportions saved to cluster_abundances_csf_annotated.csv\n")

# Visualize
p1 <- plotCellTypeProps(
  clusters = CSF$annotate,
  sample = CSF$samplename
) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  ggtitle("CSF: Cell Type Proportions")

print(p1)

## Statistical Testing for CSF
cat("\nPerforming statistical tests for CSF...\n")

df_csf <- read.csv("cluster_abundances_csf_annotated.csv")

# Reshape data for statistical testing
df_transposed_csf <- df_csf %>%
  pivot_wider(
    names_from = Counts.clusters,
    values_from = Proportions.Freq,
    id_cols = Counts.sample
  )

# Add group labels
df_transposed_csf$Group <- ifelse(
  grepl("^control_", df_transposed_csf$Counts.sample),
  "control",
  "GAD"
)

# Prepare data for volcano plot
df_processed_csf <- df_transposed_csf %>%
  select(-Counts.sample) %>%
  select(Group, everything())

# Save processed data
write.csv(
  df_processed_csf,
  "GAD_vs_control_csf_processed.csv",
  row.names = FALSE
)

## CSF Volcano Plot
cat("Generating CSF volcano plot...\n")

GAD <- df_processed_csf[df_processed_csf$Group == "GAD", -1]
control <- df_processed_csf[df_processed_csf$Group == "control", -1]

results_log2FC_csf <- data.frame(
  Cluster = colnames(GAD),
  Log2FC = NA,
  Log10pvalue = NA
)

for (i in 1:ncol(GAD)) {
  x <- GAD[, i]
  y <- control[, i]

  # Calculate fold change
  Log2FC <- log2(median(x) / median(y))

  # Calculate p-value
  pvalue <- wilcox.test(x, y)$p.value
  Log10pvalue <- -log10(pvalue)

  results_log2FC_csf[i, 2] <- Log2FC
  results_log2FC_csf[i, 3] <- Log10pvalue
}

# Remove NA and infinite values
results_log2FC_csf <- results_log2FC_csf[
  complete.cases(results_log2FC_csf$Log2FC, results_log2FC_csf$Log10pvalue) &
    !is.infinite(results_log2FC_csf$Log2FC) &
    !is.infinite(results_log2FC_csf$Log10pvalue),
]

# Save results
write.csv(results_log2FC_csf, "GAD_vs_control_CSF_Log2FC.csv", row.names = FALSE)
cat("CSF results saved to GAD_vs_control_CSF_Log2FC.csv\n")

# Create volcano plot
p <- ggplot(
  data = results_log2FC_csf,
  aes(x = Log2FC, y = Log10pvalue, label = Cluster)
) +
  geom_point(shape = 16, size = 4, aes(color = Cluster, fill = Cluster)) +
  theme_minimal() +
  geom_vline(xintercept = c(-0.485, 0.485), col = "red", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "blue", linetype = "dashed") +
  theme(
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  ) +
  geom_text_repel(aes(label = Cluster), max.overlaps = 20) +
  labs(
    x = "Log2 fold change",
    y = "-Log10 p value",
    title = "CSF: GAD vs Control"
  )

print(p)

# ---- 5.2. PBMC Compartment Analysis ----
cat("\n=== PBMC COMPARTMENT ANALYSIS ===\n")

## Cell Type Proportions in PBMC
props_pbmc <- getTransformedProps(
  clusters = PBMC$annotate,
  sample = PBMC$samplename
)

# Save proportions
write.csv(
  props_pbmc,
  "cluster_abundances_pbmc_annotated.csv",
  quote = FALSE,
  row.names = FALSE
)
cat("PBMC proportions saved to cluster_abundances_pbmc_annotated.csv\n")

# Visualize
p1 <- plotCellTypeProps(
  clusters = PBMC$annotate,
  sample = PBMC$samplename
) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  ggtitle("PBMC: Cell Type Proportions")

print(p1)

## Statistical Testing for PBMC
cat("\nPerforming statistical tests for PBMC...\n")

df_pbmc <- read.csv("cluster_abundances_pbmc_annotated.csv")

# Reshape data
df_transposed_pbmc <- df_pbmc %>%
  pivot_wider(
    names_from = Counts.clusters,
    values_from = Proportions.Freq,
    id_cols = Counts.sample
  )

# Add group labels
df_transposed_pbmc$Group <- ifelse(
  grepl("^control_", df_transposed_pbmc$Counts.sample),
  "control",
  "GAD"
)

# Prepare data
df_processed_pbmc <- df_transposed_pbmc %>%
  select(-Counts.sample) %>%
  select(Group, everything())

# Save processed data
write.csv(
  df_processed_pbmc,
  "GAD_vs_control_pbmc_processed.csv",
  row.names = FALSE
)

## PBMC Volcano Plot
cat("Generating PBMC volcano plot...\n")

GAD_pbmc <- df_processed_pbmc[df_processed_pbmc$Group == "GAD", -1]
control_pbmc <- df_processed_pbmc[df_processed_pbmc$Group == "control", -1]

results_log2FC_pbmc <- data.frame(
  Cluster = colnames(GAD_pbmc),
  Log2FC = NA,
  Log10pvalue = NA
)

for (i in 1:ncol(GAD_pbmc)) {
  x <- GAD_pbmc[, i]
  y <- control_pbmc[, i]

  Log2FC <- log2(median(x) / median(y))
  pvalue <- wilcox.test(x, y)$p.value
  Log10pvalue <- -log10(pvalue)

  results_log2FC_pbmc[i, 2] <- Log2FC
  results_log2FC_pbmc[i, 3] <- Log10pvalue
}

# Clean results
results_log2FC_pbmc <- results_log2FC_pbmc[
  complete.cases(results_log2FC_pbmc$Log2FC, results_log2FC_pbmc$Log10pvalue) &
    !is.infinite(results_log2FC_pbmc$Log2FC) &
    !is.infinite(results_log2FC_pbmc$Log10pvalue),
]

# Save results
write.csv(results_log2FC_pbmc, "GAD_vs_control_PBMC_Log2FC.csv", row.names = FALSE)
cat("PBMC results saved to GAD_vs_control_PBMC_Log2FC.csv\n")

# Create volcano plot
p <- ggplot(
  data = results_log2FC_pbmc,
  aes(x = Log2FC, y = Log10pvalue, label = Cluster)
) +
  geom_point(shape = 16, size = 4, aes(color = Cluster, fill = Cluster)) +
  theme_minimal() +
  geom_vline(xintercept = c(-0.485, 0.485), col = "red", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "blue", linetype = "dashed") +
  theme(
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  ) +
  geom_text_repel(aes(label = Cluster), max.overlaps = 20) +
  labs(
    x = "Log2 fold change",
    y = "-Log10 p value",
    title = "PBMC: GAD vs Control"
  )

print(p)

################################################################################
# 6. DIFFERENTIAL GENE EXPRESSION
################################################################################

cat("\n=== DIFFERENTIAL GENE EXPRESSION ANALYSIS ===\n")

# Set identity to disease
Idents(seurat_object_annotated) <- seurat_object_annotated$Disease

# Find DEGs
cat("Finding differentially expressed genes...\n")
dge_control_gad <- FindMarkers(
  seurat_object_annotated,
  ident.1 = "GAD",
  ident.2 = "control",
  test.use = "wilcox",
  logfc.threshold = 0.2,
  assay = "SCT"
)

# Save results
write.csv(
  dge_control_gad,
  "DEGs_GAD_vs_control_all_clusters.csv",
  row.names = TRUE
)
cat("DEG results saved to DEGs_GAD_vs_control_all_clusters.csv\n")

## Volcano Plot: All Clusters
cat("Generating volcano plot for all clusters...\n")

df_dge <- read.csv("DEGs_GAD_vs_control_all_clusters.csv", row.names = 1)
df_dge$gene_symbol <- rownames(df_dge)

vp <- EnhancedVolcano(
  df_dge,
  x = "avg_log2FC",
  y = "p_val_adj",
  lab = df_dge$gene_symbol,
  pCutoff = 1e-20,
  FCcutoff = 1.5,
  legendLabels = c(
    "Not significant",
    "Log2FC",
    "p-value",
    "p-value and Log2FC"
  ),
  legendPosition = "right",
  pointSize = 4.0,
  labSize = 5.0,
  labCol = "black",
  labFace = "bold",
  drawConnectors = TRUE,
  legendLabSize = 16,
  legendIconSize = 4.0,
  title = "GAD versus Control (All Clusters)"
)

print(vp)

################################################################################
# 7. PSEUDOBULK VOLCANO PLOT
################################################################################

cat("\n=== PSEUDOBULK DIFFERENTIAL EXPRESSION ANALYSIS ===\n")

# ---- Configuration ----
LOG2FC_THRESHOLD <- 0.585  # ~1.5-fold change
PADJ_THRESHOLD <- 0.05
N_TOP_GENES <- 10  # Number of top genes to label
INPUT_FILE <- "pseudobulk_deseq2_GAD_CSF_expanded_vs_nonexpanded.csv"

# ---- Load differential expression results ----
cat(sprintf("Loading results from %s...\n", INPUT_FILE))
de_results <- read.csv(INPUT_FILE, row.names = 1) %>%
  rownames_to_column(var = "gene") %>%
  as_tibble()

cat(sprintf("Loaded %d genes\n", nrow(de_results)))

# ---- Define color scheme for significance ----
color_scheme <- case_when(
  de_results$log2FoldChange <= -LOG2FC_THRESHOLD & de_results$padj <= PADJ_THRESHOLD ~ "dodgerblue3",
  de_results$log2FoldChange >= LOG2FC_THRESHOLD & de_results$padj <= PADJ_THRESHOLD ~ "chartreuse3",
  TRUE ~ "grey"
)

# Replace NA values with grey
color_scheme[is.na(color_scheme)] <- "grey"

# Create named vector for legend
names(color_scheme)[color_scheme == "chartreuse3"] <- "Up-regulated"
names(color_scheme)[color_scheme == "dodgerblue3"] <- "Down-regulated"
names(color_scheme)[color_scheme == "grey"] <- "Not significant"

# ---- Identify top differentially expressed genes ----
cat("Identifying top differentially expressed genes...\n")

# Top up-regulated genes
top_upregulated <- de_results %>%
  filter(log2FoldChange >= LOG2FC_THRESHOLD, padj <= PADJ_THRESHOLD) %>%
  arrange(padj, desc(log2FoldChange)) %>%
  head(N_TOP_GENES)

# Top down-regulated genes
top_downregulated <- de_results %>%
  filter(log2FoldChange <= -LOG2FC_THRESHOLD, padj <= PADJ_THRESHOLD) %>%
  arrange(padj, log2FoldChange) %>%
  head(N_TOP_GENES)

# Combine top genes for labeling
genes_to_label <- c(top_upregulated$gene, top_downregulated$gene)

cat(sprintf("\nIdentified %d up-regulated and %d down-regulated genes\n",
            nrow(top_upregulated), nrow(top_downregulated)))
cat("Top labeled genes:", paste(genes_to_label, collapse = ", "), "\n\n")

# ---- Generate volcano plot ----
cat("Generating enhanced volcano plot...\n")

volcano_plot <- EnhancedVolcano(
  de_results,
  lab = de_results$gene,
  x = "log2FoldChange",
  y = "padj",

  # Significance thresholds
  pCutoff = PADJ_THRESHOLD,
  FCcutoff = LOG2FC_THRESHOLD,

  # Point styling
  pointSize = 3,
  colAlpha = 0.6,
  colCustom = color_scheme,
  shape = 16,

  # Axis labels
  xlab = bquote(~Log[2]~ "fold change"),
  ylab = bquote(~-Log[10]~ "adjusted p-value"),

  # Gene labels
  selectLab = genes_to_label,
  boxedLabels = TRUE,
  labSize = 3.5,
  labFace = "bold",
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  lengthConnectors = unit(0.01, "npc"),
  arrowheads = FALSE,
  maxoverlapsConnectors = Inf,
  max.overlaps = Inf,

  # Plot annotations
  title = "Differential Expression: Expanded vs Non-Expanded Clones",
  subtitle = "GAD Patient CSF T cells (Pseudobulk Analysis)",
  caption = bquote(~Log[2]~ "FC cutoff: ±" * .(LOG2FC_THRESHOLD) * 
                   "; Adjusted p-value cutoff: " * .(PADJ_THRESHOLD)),

  # Legend
  legendPosition = "right",
  legendLabSize = 12,
  legendIconSize = 4.0
)

print(volcano_plot)

# ---- Summary statistics ----
summary_stats <- de_results %>%
  mutate(
    regulation = case_when(
      log2FoldChange >= LOG2FC_THRESHOLD & padj <= PADJ_THRESHOLD ~ "Up-regulated",
      log2FoldChange <= -LOG2FC_THRESHOLD & padj <= PADJ_THRESHOLD ~ "Down-regulated",
      TRUE ~ "Not significant"
    )
  ) %>%
  count(regulation)

cat("\n=== Differential Expression Summary ===\n")
print(summary_stats)


################################################################################
# 8. SESSION INFORMATION
################################################################################

cat("\n=== SESSION INFORMATION ===\n")
sessionInfo()

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("All results have been saved to the working directory.\n")
