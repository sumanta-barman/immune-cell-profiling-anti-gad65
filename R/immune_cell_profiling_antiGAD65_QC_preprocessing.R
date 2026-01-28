################################################################################
# Immune cell profiling anti-GAD65 – QC and preprocessing
# Author: Sumanta Barman
# Date: 2026-01-26
# Description: Single-cell RNA-seq QC and preprocessing pipeline including
#              quality control, doublet removal, normalization, and batch correction
################################################################################

# Load required libraries --------------------------------------------------
cat("Loading required libraries...\n")
suppressPackageStartupMessages({
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
})

# Set system environment
Sys.setenv(lang = "en_US")
options(Seurat.object.assay.version = "v5")
theme_set(theme_classic())

cat("Libraries loaded successfully.\n\n")

# Load raw data and initial QC ---------------------------------------------
cat("=== SECTION 1: Loading raw data ===\n")

# Define data directory
data_dir <- "/Data/"

# List all sample directories
dirs <- list.dirs(path = data_dir, recursive = FALSE, full.names = FALSE)
cat("Found", length(dirs), "sample directories:\n")
print(dirs)

# Initialize containers for Seurat objects
seurat_objects <- list()
seurat_object_names <- character(0)

# Load 10X data for each sample
for (x in dirs) {
  sample_name <- gsub("_filtered_feature_bc_matrix", "", x)

  cat("Loading", sample_name, "...\n")
  counts <- Read10X(data.dir = file.path(data_dir, x))
  seurat_obj <- CreateSeuratObject(counts = counts)

  seurat_objects[[sample_name]] <- seurat_obj
  seurat_object_names <- c(seurat_object_names, sample_name)
}

cat("\nLoaded", length(seurat_objects), "samples\n\n")

# Define QC analysis function ----------------------------------------------
cat("=== SECTION 2: Initial QC analysis ===\n")

analyze_seurat_object <- function(seurat_obj, filename_prefix) {

  # Add sample metadata
  seurat_obj$sample <- filename_prefix

  # Calculate QC metrics
  seurat_obj$mitoPercent <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj$riboPercent <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")
  seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)

  # Generate violin plots
  p_vln <- VlnPlot(
    seurat_obj,
    features = c("nFeature_RNA", "nCount_RNA", "mitoPercent", "riboPercent"),
    ncol = 4
  )
  ggsave(
    filename = paste0(filename_prefix, "_QC_ViolinPlot.pdf"),
    plot = p_vln,
    device = "pdf"
  )

  # Generate scatter plots
  scatter_plot1 <- FeatureScatter(
    seurat_obj,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA"
  ) + geom_smooth(method = "lm")

  scatter_plot2 <- FeatureScatter(
    seurat_obj,
    feature1 = "nCount_RNA",
    feature2 = "mitoPercent"
  ) + geom_smooth(method = "lm")

  scatter_plot3 <- FeatureScatter(
    seurat_obj,
    feature1 = "nCount_RNA",
    feature2 = "riboPercent"
  ) + geom_smooth(method = "lm")

  scatter_combined <- CombinePlots(
    plots = list(scatter_plot1, scatter_plot2, scatter_plot3),
    ncol = 1
  )
  ggsave(
    filename = paste0(filename_prefix, "_QC_ScatterPlot.pdf"),
    plot = scatter_combined,
    device = "pdf"
  )

  # Export metadata
  metadata <- as.data.frame(seurat_obj@meta.data)
  write.csv(
    metadata,
    file = paste0("QC_metadata_", filename_prefix, ".csv"),
    row.names = TRUE
  )

  # Calculate outlier thresholds using scuttle
  qc_mito <- isOutlier(seurat_obj$mitoPercent, log = FALSE, type = "higher")
  qc_ribo <- isOutlier(seurat_obj$riboPercent, log = FALSE, type = "higher")
  qc_nFeature_higher <- isOutlier(seurat_obj$nFeature_RNA, log = TRUE, type = "higher")
  qc_nFeature_lower <- isOutlier(seurat_obj$nFeature_RNA, log = TRUE, type = "lower")

  qc_info <- data.frame(
    mitoPercent = attr(qc_mito, "thresholds"),
    riboPercent = attr(qc_ribo, "thresholds"),
    nFeature_RNA_higher = attr(qc_nFeature_higher, "thresholds"),
    nFeature_RNA_lower = attr(qc_nFeature_lower, "thresholds")
  )

  write.csv(
    qc_info,
    file = paste0("QC_info_", filename_prefix, ".csv"),
    row.names = TRUE
  )

  # Add cell IDs to metadata
  metadata$cells <- rownames(metadata)
  seurat_obj@meta.data <- metadata

  # Generate density plots
  p1 <- metadata %>%
    ggplot(aes(x = sample, fill = sample)) +
    geom_bar() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    ggtitle("Number of Cells")

  p2 <- metadata %>%
    ggplot(aes(color = sample, x = nCount_RNA, fill = sample)) +
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    ylab("Cell density") +
    geom_vline(xintercept = 500) +
    ggtitle("nCount_RNA (UMI) > 500")

  p3 <- metadata %>%
    ggplot(aes(color = sample, x = nFeature_RNA, fill = sample)) +
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    geom_vline(xintercept = 300) +
    ggtitle("nFeature_RNA (genes) > 300")

  p4 <- metadata %>%
    ggplot(aes(x = nCount_RNA, y = nFeature_RNA, color = mitoPercent)) +
    geom_point() +
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method = "lm") +
    scale_x_log10() +
    scale_y_log10() +
    geom_vline(xintercept = 500) +
    geom_hline(yintercept = 250) +
    facet_wrap(~sample) +
    ggtitle("nUMI vs nGenes")

  p5 <- metadata %>%
    ggplot(aes(color = sample, x = mitoPercent, fill = sample)) +
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    geom_vline(xintercept = 10) +
    ggtitle("Mitochondrial genes < 10%")

  p6 <- metadata %>%
    ggplot(aes(x = log10GenesPerUMI, color = sample, fill = sample)) +
    geom_density(alpha = 0.2) +
    geom_vline(xintercept = 0.8) +
    ggtitle("Complexity (log10GenesPerUMI)")

  p7 <- metadata %>%
    ggplot(aes(color = sample, x = riboPercent, fill = sample)) +
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    geom_vline(xintercept = 60) +
    ggtitle("Ribosomal genes < 60%")

  p8 <- metadata %>%
    ggplot(aes(color = sample, x = nFeature_RNA, fill = sample)) +
    geom_density(alpha = 0.2) +
    scale_x_log10() +
    geom_vline(xintercept = 3000) +
    ggtitle("nFeature_RNA < 3000")

  # Optional: save combined panels
  plot_combined1 <- p2 + p3 + p8 + plot_layout(ncol = 1)
  plot_combined2 <- p5 + p7 + plot_layout(ncol = 1)
  plot_combined3 <- p1 + p4 + p6 + plot_layout(ncol = 1)

  return(seurat_obj)
}

# Run QC analysis for all samples ------------------------------------------
for (sample_name in seurat_object_names) {
  cat("Analyzing", sample_name, "...\n")
  seurat_objects[[sample_name]] <- analyze_seurat_object(
    seurat_objects[[sample_name]],
    sample_name
  )
}

cat("\nQC analysis completed for all samples\n\n")

# Display summary of all Seurat objects
print(seurat_objects)

# Sample-specific QC filtering ---------------------------------------------
cat("\n=== SECTION 3: Sample-specific QC filtering ===\n")

# CSF Control samples
cat("Filtering CSF Control samples...\n")
control_CSF_8360_filtered <- subset(
  seurat_objects$Control_CSF_8360,
  subset = nFeature_RNA > 200 & nFeature_RNA < 3000 &
    mitoPercent < 8 & nCount_RNA > 500 & log10GenesPerUMI > 0.70
)

control_CSF_8361_filtered <- subset(
  seurat_objects$Control_CSF_8361,
  subset = nFeature_RNA > 200 & nFeature_RNA < 3000 &
    mitoPercent < 8 & nCount_RNA > 500 & log10GenesPerUMI > 0.80
)

control_CSF_MPST83775_filtered <- subset(
  seurat_objects$Control_CSF_MPST83775,
  subset = nFeature_RNA > 300 & nFeature_RNA < 3000 &
    mitoPercent < 7 & nCount_RNA > 500 & log10GenesPerUMI > 0.80
)

control_CSF_MPST95809_filtered <- subset(
  seurat_objects$Control_CSF_MPST95809,
  subset = nFeature_RNA > 250 & nFeature_RNA < 3000 &
    mitoPercent < 8 & nCount_RNA > 500 & log10GenesPerUMI > 0.80
)

control_CSF_R21_filtered <- subset(
  seurat_objects$Control_CSF_R21,
  subset = nFeature_RNA > 300 & nFeature_RNA < 3000 &
    mitoPercent < 5 & nCount_RNA > 500 & log10GenesPerUMI > 0.80
)

control_CSF_R31_filtered <- subset(
  seurat_objects$Control_CSF_R31,
  subset = nFeature_RNA > 300 & nFeature_RNA < 3500 &
    mitoPercent < 5 & nCount_RNA > 500 & log10GenesPerUMI > 0.80
)

control_CSF_R32_filtered <- subset(
  seurat_objects$Control_CSF_R32,
  subset = nFeature_RNA > 300 & nFeature_RNA < 3500 &
    mitoPercent < 6 & nCount_RNA > 500 & log10GenesPerUMI > 0.80
)

control_CSF_R34_filtered <- subset(
  seurat_objects$Control_CSF_R34,
  subset = nFeature_RNA > 500 & nFeature_RNA < 3500 &
    mitoPercent < 6 & nCount_RNA > 1000 & log10GenesPerUMI > 0.80
)

# GAD CSF samples
cat("Filtering GAD CSF samples...\n")
GAD_CSF_9910_filtered <- subset(
  seurat_objects$GAD_CSF_9910,
  subset = nFeature_RNA > 400 & nFeature_RNA < 4500 &
    mitoPercent < 9 & nCount_RNA > 500 & log10GenesPerUMI > 0.80
)

GAD_CSF_9961_filtered <- subset(
  seurat_objects$GAD_CSF_9961,
  subset = nFeature_RNA > 250 & nFeature_RNA < 3500 &
    mitoPercent < 12 & nCount_RNA > 400 & log10GenesPerUMI > 0.80
)

GAD_CSF_A15_filtered <- subset(
  seurat_objects$GAD_CSF_A15,
  subset = nFeature_RNA > 120 & nFeature_RNA < 750 &
    mitoPercent < 5 & nCount_RNA > 300 & log10GenesPerUMI > 0.80
)

GAD_CSF_A18_filtered <- subset(
  seurat_objects$GAD_CSF_A18,
  subset = nFeature_RNA > 100 & nFeature_RNA < 700 &
    mitoPercent < 6 & nCount_RNA > 200 & log10GenesPerUMI > 0.80
)

GAD_CSF_A24_filtered <- subset(
  seurat_objects$GAD_CSF_A24,
  subset = nFeature_RNA > 300 & nFeature_RNA < 3000 &
    mitoPercent < 6 & nCount_RNA > 500 & log10GenesPerUMI > 0.80
)

GAD_CSF_R3_filtered <- subset(
  seurat_objects$GAD_CSF_R3,
  subset = nFeature_RNA > 250 & nFeature_RNA < 3000 &
    mitoPercent < 7 & nCount_RNA > 500 & log10GenesPerUMI > 0.80
)

GAD_CSF_R14_filtered <- subset(
  seurat_objects$GAD_CSF_R14,
  subset = nFeature_RNA > 300 & nFeature_RNA < 3000 &
    mitoPercent < 8 & nCount_RNA > 500 & log10GenesPerUMI > 0.80
)

# PBMC Control samples
cat("Filtering PBMC Control samples...\n")
control_PBMC_8360_filtered <- subset(
  seurat_objects$Control_PBMC_8360,
  subset = nFeature_RNA > 300 & nFeature_RNA < 3000 &
    mitoPercent < 12 & nCount_RNA > 500 & log10GenesPerUMI > 0.80
)

control_PBMC_8361_filtered <- subset(
  seurat_objects$Control_PBMC_8361,
  subset = nFeature_RNA > 300 & nFeature_RNA < 3000 &
    mitoPercent < 10 & nCount_RNA > 500 & log10GenesPerUMI > 0.80
)

control_PBMC_MPST83775_filtered <- subset(
  seurat_objects$Control_PBMC_MPST83775,
  subset = nFeature_RNA > 400 & nFeature_RNA < 3000 &
    mitoPercent < 7 & nCount_RNA > 500 & log10GenesPerUMI > 0.80
)

control_PBMC_MPST95809_filtered <- subset(
  seurat_objects$Control_PBMC_MPST95809,
  subset = nFeature_RNA > 300 & nFeature_RNA < 3000 &
    mitoPercent < 10 & nCount_RNA > 500 & log10GenesPerUMI > 0.80
)

control_PBMC_R21_filtered <- subset(
  seurat_objects$Control_PBMC_R21,
  subset = nFeature_RNA > 300 & nFeature_RNA < 3000 &
    mitoPercent < 6 & nCount_RNA > 500 & log10GenesPerUMI > 0.80
)

control_PBMC_R31_filtered <- subset(
  seurat_objects$Control_PBMC_R31,
  subset = nFeature_RNA > 300 & nFeature_RNA < 3000 &
    mitoPercent < 7 & nCount_RNA > 500 & log10GenesPerUMI > 0.80
)

control_PBMC_R32_filtered <- subset(
  seurat_objects$Control_PBMC_R32,
  subset = nFeature_RNA > 400 & nFeature_RNA < 3500 &
    mitoPercent < 6 & nCount_RNA > 500 & log10GenesPerUMI > 0.80
)

control_PBMC_R34_filtered <- subset(
  seurat_objects$Control_PBMC_R34,
  subset = nFeature_RNA > 300 & nFeature_RNA < 3500 &
    mitoPercent < 5 & nCount_RNA > 500 & log10GenesPerUMI > 0.80
)

# GAD PBMC samples
cat("Filtering GAD PBMC samples...\n")
GAD_PBMC_9910_filtered <- subset(
  seurat_objects$GAD_PBMC_9910,
  subset = nFeature_RNA > 300 & nFeature_RNA < 3000 &
    mitoPercent < 20 & nCount_RNA > 500 & log10GenesPerUMI > 0.80
)

GAD_PBMC_9961_filtered <- subset(
  seurat_objects$GAD_PBMC_9961,
  subset = nFeature_RNA > 300 & nFeature_RNA < 3000 &
    mitoPercent < 20 & nCount_RNA > 500 & log10GenesPerUMI > 0.80
)

GAD_PBMC_KKS_filtered <- subset(
  seurat_objects$GAD_PBMC_KKS,
  subset = nFeature_RNA > 300 & nFeature_RNA < 4000 &
    mitoPercent < 9 & nCount_RNA > 500 & log10GenesPerUMI > 0.80
)

GAD_PBMC_R3_filtered <- subset(
  seurat_objects$GAD_PBMC_R3,
  subset = nFeature_RNA > 300 & nFeature_RNA < 3000 &
    mitoPercent < 8 & nCount_RNA > 500 & log10GenesPerUMI > 0.80
)

GAD_PBMC_R14_filtered <- subset(
  seurat_objects$GAD_PBMC_R14,
  subset = nFeature_RNA > 300 & nFeature_RNA < 3000 &
    mitoPercent < 10 & nCount_RNA > 500 & log10GenesPerUMI > 0.80
)

# Collect all filtered objects
filtered_object_names <- c(
  "control_CSF_8360_filtered", "control_CSF_8361_filtered",
  "control_CSF_MPST83775_filtered", "control_CSF_MPST95809_filtered",
  "control_CSF_R21_filtered", "control_CSF_R31_filtered",
  "control_CSF_R32_filtered", "control_CSF_R34_filtered",
  "GAD_CSF_9910_filtered", "GAD_CSF_9961_filtered",
  "GAD_CSF_A15_filtered", "GAD_CSF_A18_filtered",
  "GAD_CSF_A24_filtered", "GAD_CSF_R3_filtered", "GAD_CSF_R14_filtered",
  "control_PBMC_8360_filtered", "control_PBMC_8361_filtered",
  "control_PBMC_MPST83775_filtered", "control_PBMC_MPST95809_filtered",
  "control_PBMC_R21_filtered", "control_PBMC_R31_filtered",
  "control_PBMC_R32_filtered", "control_PBMC_R34_filtered",
  "GAD_PBMC_9910_filtered", "GAD_PBMC_9961_filtered",
  "GAD_PBMC_KKS_filtered", "GAD_PBMC_R3_filtered", "GAD_PBMC_R14_filtered"
)

cat("\nFiltered", length(filtered_object_names), "samples\n\n")

# Post-filtering QC visualization ------------------------------------------
cat("=== SECTION 4: Post-filtering QC visualization ===\n")

generate_plots <- function(data, filename_prefix) {
  # Violin plots
  p_vln <- VlnPlot(
    data,
    features = c("nFeature_RNA", "nCount_RNA", "mitoPercent", "riboPercent"),
    ncol = 4
  )
  ggsave(
    filename = paste0(filename_prefix, "_ViolinPlot.pdf"),
    plot = p_vln,
    device = "pdf",
    dpi = 700
  )

  # Scatter plots
  scatter_plot1 <- FeatureScatter(
    data,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA"
  ) + geom_smooth(method = "lm")

  scatter_plot2 <- FeatureScatter(
    data,
    feature1 = "nCount_RNA",
    feature2 = "mitoPercent"
  ) + geom_smooth(method = "lm")

  scatter_plot3 <- FeatureScatter(
    data,
    feature1 = "nCount_RNA",
    feature2 = "riboPercent"
  ) + geom_smooth(method = "lm")

  scatter_combined <- CombinePlots(
    plots = list(scatter_plot1, scatter_plot2, scatter_plot3),
    ncol = 1
  )
  ggsave(
    filename = paste0(filename_prefix, "_ScatterPlot.pdf"),
    plot = scatter_combined,
    device = "pdf",
    dpi = 700
  )
}

# Generate plots for all filtered samples
for (obj_name in filtered_object_names) {
  cat("Generating post-QC plots for", obj_name, "\n")
  generate_plots(get(obj_name), paste0("PostQC_", obj_name))
}

cat("\nPost-filtering QC visualization completed\n\n")

# Doublet detection and removal --------------------------------------------
cat("=== SECTION 5: Doublet detection and removal ===\n")

# Get all filtered objects
filtered_objects <- mget(filtered_object_names)

# Merge all Seurat objects
cat("Merging all filtered samples...\n")
combined <- merge(
  x = filtered_objects[[1]],
  y = filtered_objects[2:length(filtered_objects)],
  add.cell.ids = filtered_object_names,
  project = "Immune_cell_profiling_anti_GAD65"
)

cat("Combined object created:\n")
print(combined)

saveRDS(combined, file = "seurat_object_combined.rds")
cat("Saved: seurat_object_combined.rds\n\n")

# Display sample distribution
cat("Sample distribution:\n")
print(table(combined$sample))
cat("\n")

# Split combined object by sample
cat("Splitting by sample for DoubletFinder...\n")
combined.split <- SplitObject(combined, split.by = "sample")

# Run DoubletFinder per sample
for (i in seq_along(combined.split)) {
  sample_name <- names(combined.split)[i]
  cat("\n--- Processing sample", i, "of", length(combined.split), ":", sample_name, "---\n")

  seu_sample <- combined.split[[i]]

  # SCTransform normalization
  cat("  Running SCTransform...\n")
  seu_sample <- SCTransform(seu_sample, verbose = FALSE)

  # PCA
  cat("  Running PCA...\n")
  seu_sample <- RunPCA(seu_sample, npcs = 50, verbose = FALSE)

  # Determine optimal number of PCs
  stdv <- seu_sample[["pca"]]@stdev
  pct <- (stdv / sum(stdv)) * 100
  cumulative <- cumsum(pct)

  co1 <- which(cumulative > 90 & pct < 5)[1]
  co2 <- sort(
    which((pct[1:(length(pct) - 1)] - pct[2:length(pct)]) > 0.1),
    decreasing = TRUE
  )[1] + 1

  min_pc <- min(co1, co2, na.rm = TRUE)
  if (is.na(min_pc) || is.infinite(min_pc)) {
    min_pc <- 20
  }

  cat("  Using", min_pc, "PCs\n")

  # UMAP and clustering
  cat("  Running UMAP and clustering...\n")
  seu_sample <- RunUMAP(seu_sample, dims = 1:min_pc, verbose = FALSE)
  seu_sample <- FindNeighbors(seu_sample, dims = 1:min_pc, verbose = FALSE)
  seu_sample <- FindClusters(seu_sample, resolution = 0.8, verbose = FALSE)

  # DoubletFinder parameter sweep
  cat("  Running DoubletFinder parameter sweep...\n")
  sweep.list <- paramSweep(seu_sample, PCs = 1:min_pc, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)

  # Find optimal pK
  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric), ]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(as.character(optimal.pk))

  cat("  Optimal pK:", optimal.pk, "\n")

  # Estimate doublet rate
  annotations <- seu_sample@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp.poi <- round(optimal.pk * nrow(seu_sample@meta.data))
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))

  cat("  Expected doublets:", nExp.poi.adj, "\n")

  # Run DoubletFinder
  cat("  Running DoubletFinder...\n")
  seu_sample <- doubletFinder(
    seu = seu_sample,
    PCs = 1:min_pc,
    pK = optimal.pk,
    nExp = nExp.poi.adj,
    sct = TRUE
  )

  # Rename doublet classification column
  metadata <- seu_sample@meta.data
  doublet_col <- grep("DF.classifications", colnames(metadata), value = TRUE)
  colnames(metadata)[colnames(metadata) == doublet_col] <- "doublet_finder"
  seu_sample@meta.data <- metadata

  # Subset singlets only
  seu_singlets <- subset(seu_sample, subset = doublet_finder == "Singlet")

  cat("  Retained", ncol(seu_singlets), "singlets out of", ncol(seu_sample), "cells\n")

  combined.split[[i]] <- seu_singlets
}

# Merge all singlet objects
cat("\nMerging all singlet objects...\n")
combined_singlets <- Reduce(
  f = function(x, y) merge(x, y, project = "GAD_control_PBMC_CSF"),
  x = combined.split
)

cat("Combined singlets object created:\n")
print(combined_singlets)

saveRDS(combined_singlets, file = "seurat_object_combined_singlets.rds")
cat("Saved: seurat_object_combined_singlets.rds\n\n")

# Add metadata columns -----------------------------------------------------
cat("=== SECTION 6: Adding metadata columns ===\n")

# Extract sample names
samples <- combined_singlets@meta.data$sample

# Split sample names into components
split_samples <- strsplit(as.character(samples), "_")

# Create new metadata columns
new_columns <- data.frame(
  Disease = sapply(split_samples, function(x) x[1]),
  Sampletype = sapply(split_samples, function(x) x[2]),
  Patient = sapply(split_samples, function(x) x[3]),
  Center = sapply(split_samples, function(x) ifelse(length(x) >= 4, x[4], NA)),
  row.names = rownames(combined_singlets@meta.data)
)

# Add new columns to metadata
combined_singlets@meta.data <- cbind(combined_singlets@meta.data, new_columns)

# Create combined sample name
combined_singlets$samplename <- paste(
  combined_singlets$Disease,
  combined_singlets$Sampletype,
  combined_singlets$Patient,
  sep = "_"
)

cat("Added metadata columns: Disease, Sampletype, Patient, Center, samplename\n")
cat("\nFirst few rows of metadata:\n")
print(head(combined_singlets@meta.data))
cat("\n")

# Normalization and dimensionality reduction -------------------------------
cat("=== SECTION 7: Normalization and dimensionality reduction ===\n")

# SCTransform normalization
cat("Running SCTransform normalization...\n")
combined_singlets <- SCTransform(
  object = combined_singlets,
  vst.flavor = "v2",
  method = "glmGamPoi",
  verbose = FALSE
)

# Run PCA
cat("Running PCA...\n")
combined_singlets <- RunPCA(combined_singlets, npcs = 50, verbose = FALSE)

# Determine optimal number of PCs
cat("\nDetermining optimal number of PCs...\n")

# Calculate PC statistics
pct <- combined_singlets[["pca"]]@stdev / sum(combined_singlets[["pca"]]@stdev) * 100
cumu <- cumsum(pct)

# Threshold 1: cumulative > 90% and individual PC < 5%
co1 <- which(cumu > 90 & pct < 5)[1]
cat("PC threshold (90% cumulative, <5% individual):", co1, "\n")

# Threshold 2: PC variance change < 0.1%
co2 <- sort(
  which((pct[1:(length(pct) - 1)] - pct[2:length(pct)]) > 0.1),
  decreasing = TRUE
)[1] + 1
cat("PC threshold (variance change <0.1%):", co2, "\n")

# Use minimum of two thresholds
pcs <- min(co1, co2, na.rm = TRUE)
if (is.na(pcs) || is.infinite(pcs)) {
  pcs <- 20
}
cat("Using", pcs, "PCs for downstream analysis\n\n")

# Pre-integration clustering -----------------------------------------------
cat("=== SECTION 8: Pre-integration clustering ===\n")

# Clustering before batch correction
cat("Running clustering before integration...\n")
combined_singlets <- FindNeighbors(combined_singlets, dims = 1:pcs)
combined_singlets <- FindClusters(combined_singlets, resolution = 0.5)
combined_singlets <- RunUMAP(combined_singlets, dims = 1:pcs)

# Visualize batch effects
cat("Generating pre-integration UMAP plots...\n")
p1 <- DimPlot(combined_singlets, reduction = "umap", group.by = "sample") +
  ggtitle("Before Integration - by Sample")
ggsave("Pre_integration_by_sample.pdf", plot = p1, width = 12, height = 8)

p2 <- DimPlot(combined_singlets, reduction = "umap", group.by = "Disease") +
  ggtitle("Before Integration - by Disease")
ggsave("Pre_integration_by_disease.pdf", plot = p2, width = 10, height = 8)

p3 <- DimPlot(combined_singlets, reduction = "umap", group.by = "Center") +
  ggtitle("Before Integration - by Center")
ggsave("Pre_integration_by_center.pdf", plot = p3, width = 10, height = 8)

cat("Pre-integration plots saved\n\n")

# Batch correction using Harmony -------------------------------------------
cat("=== SECTION 9: Batch correction using Harmony ===\n")

# Integrate using Harmony
cat("Running Harmony integration...\n")
combined_singlets_integrated <- IntegrateLayers(
  object = combined_singlets,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  verbose = FALSE,
  assay = "SCT"
)

saveRDS(combined_singlets_integrated, file = "seurat_object_combined_singlets_integrated.rds")
cat("Saved: seurat_object_combined_singlets_integrated.rds\n\n")

# Post-integration analysis ------------------------------------------------
cat("=== SECTION 10: Post-integration analysis and visualization ===\n")

# Clustering on harmony embeddings
cat("Running clustering on Harmony embeddings...\n")
combined_singlets_integrated <- FindNeighbors(
  combined_singlets_integrated,
  reduction = "harmony",
  dims = 1:20
)

combined_singlets_integrated <- FindClusters(
  combined_singlets_integrated,
  resolution = 0.5,
  cluster.name = "harmony_clusters"
)

# UMAP on harmony embeddings
cat("Running UMAP on Harmony embeddings...\n")
combined_singlets_integrated <- RunUMAP(
  combined_singlets_integrated,
  reduction = "harmony",
  dims = 1:20,
  reduction.name = "umap.harmony"
)

# Visualization of integrated data -----------------------------------------
cat("\nGenerating post-integration visualizations...\n")

# Clusters with labels
p_clusters <- DimPlot(
  combined_singlets_integrated,
  reduction = "umap.harmony",
  group.by = "harmony_clusters",
  label = TRUE,
  label.size = 3,
  raster = FALSE
) + ggtitle("Harmony Clusters")
ggsave("Post_integration_clusters.pdf", plot = p_clusters, width = 10, height = 8)

# By sample
p_sample <- DimPlot(
  combined_singlets_integrated,
  reduction = "umap.harmony",
  group.by = "sample"
) + ggtitle("Integrated - by Sample")
ggsave("Post_integration_by_sample.pdf", plot = p_sample, width = 12, height = 8)

# By disease
p_disease <- DimPlot(
  combined_singlets_integrated,
  reduction = "umap.harmony",
  group.by = "Disease"
) + ggtitle("Integrated - by Disease")
ggsave("Post_integration_by_disease.pdf", plot = p_disease, width = 10, height = 8)

# By sample type
p_sampletype <- DimPlot(
  combined_singlets_integrated,
  reduction = "umap.harmony",
  group.by = "Sampletype"
) + ggtitle("Integrated - by Sample Type")
ggsave("Post_integration_by_sampletype.pdf", plot = p_sampletype, width = 10, height = 8)

# By center
p_center <- DimPlot(
  combined_singlets_integrated,
  reduction = "umap.harmony",
  group.by = "Center"
) + ggtitle("Integrated - by Center")
ggsave("Post_integration_by_center.pdf", plot = p_center, width = 10, height = 8)

cat("Post-integration plots saved\n\n")

# Session information ------------------------------------------------------
cat("=== Analysis completed successfully! ===\n\n")
cat("Session Information:\n")
print(sessionInfo())

cat("\n=== END OF SCRIPT ===\n")
