################################################################################
# Package Installation Script
# Install all required packages for immune cell profiling analysis
################################################################################

# Function to check and install packages
install_if_missing <- function(pkg, bioconductor = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (bioconductor) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg, update = FALSE)
    } else {
      install.packages(pkg)
    }
    cat(paste("Installed:", pkg, "\n"))
  } else {
    cat(paste("Already installed:", pkg, "\n"))
  }
}

# CRAN packages
cran_packages <- c(
  "tidyverse", "ggplot2", "dplyr", "tidyr", "readr",
  "patchwork", "cowplot", "gridExtra", "rafalib",
  "ggrepel", "RColorBrewer", "ggsci", "ggpubr",
  "rstatix", "pheatmap", "Matrix",
  "reticulate", "devtools", "remotes"
)

# Bioconductor packages
bioc_packages <- c(
  "Seurat", "SeuratDisk", "SeuratData", "SeuratWrappers",
  "SingleR", "celldex", "Azimuth", "Signac",
  "harmony", "sctransform", "scuttle",
  "limma", "edgeR", "speckle",
  "DoubletFinder"
)

# GitHub packages
github_packages <- list(
  "immunarch" = "immunomind/immunarch",
  "EnhancedVolcano" = "kevinblighe/EnhancedVolcano"
)

cat("Installing CRAN packages...\n")
for (pkg in cran_packages) {
  install_if_missing(pkg, bioconductor = FALSE)
}

cat("\nInstalling Bioconductor packages...\n")
for (pkg in bioc_packages) {
  install_if_missing(pkg, bioconductor = TRUE)
}

cat("\nInstalling GitHub packages...\n")
for (pkg_name in names(github_packages)) {
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    remotes::install_github(github_packages[[pkg_name]])
    cat(paste("Installed:", pkg_name, "from GitHub\n"))
  } else {
    cat(paste("Already installed:", pkg_name, "\n"))
  }
}

cat("\n=== Package installation complete! ===\n")
cat("\nVerifying key packages...\n")

# Verify installations
key_packages <- c("Seurat", "SingleR", "Azimuth", "speckle", "immunarch", "EnhancedVolcano")
for (pkg in key_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(paste("✓", pkg, "successfully installed\n"))
  } else {
    cat(paste("✗", pkg, "FAILED to install\n"))
  }
}

cat("\nAll done! You can now run the analysis.\n")


################################################################################
# Functions for Immune Cell Profiling Analysis
# Author: Sumanta Barman
# Purpose: Companion functions for scRNA-seq analysis of anti-GAD65 encephalitis
################################################################################

#' Custom flat violin plot geometry for ggplot2
#'
#' @description Creates half-violin plots for visualizing distributions
#' @export
"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}

geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}

#' GeomFlatViolin ggproto object
#' @format NULL
#' @usage NULL
#' @export
GeomFlatViolin <- ggproto(
  "GeomFlatViolin", Geom,
  setup_data = function(data, params) {
    data$width <- data$width %||%
      params$width %||% (resolution(data$x, FALSE) * 0.9)
    
    data %>%
      group_by(group) %>%
      mutate(
        ymin = min(y),
        ymax = max(y),
        xmin = x,
        xmax = x + width / 2
      )
  },
  
  draw_group = function(data, panel_scales, coord) {
    data <- transform(
      data, 
      xminv = x,
      xmaxv = x + violinwidth * (xmax - x)
    )
    
    newdata <- rbind(
      plyr::arrange(transform(data, x = xminv), y),
      plyr::arrange(transform(data, x = xmaxv), -y)
    )
    
    newdata <- rbind(newdata, newdata[1,])
    
    ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord))
  },
  
  draw_key = draw_key_polygon,
  
  default_aes = aes(
    weight = 1, 
    colour = "grey20", 
    fill = "white", 
    size = 0.5,
    alpha = NA, 
    linetype = "solid"
  ),
  
  required_aes = c("x", "y")
)

#' Create raincloud plot for cell type proportions
#'
#' @param data Data frame with proportions
#' @param group_col Column name for grouping variable
#' @param value_col Column name for values to plot
#' @param colors Named vector of colors for groups
#' @param ylim Y-axis limits
#' @return ggplot object
#' @export
plot_raincloud <- function(data, group_col, value_col, colors, ylim = c(-0.2, 1.6)) {
  require(ggplot2)
  require(cowplot)
  
  plotData <- data %>%
    mutate(
      value = as.numeric(.data[[value_col]]),
      Group = as.factor(.data[[group_col]])
    ) %>%
    filter(!is.na(value))
  
  ggplot(plotData, aes(x = Group, y = value, fill = Group, color = Group, alpha = 0.8)) +
    theme_cowplot() +
    ylim(ylim) +
    scale_shape_identity() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 12),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 10),
      axis.text.x = element_text(angle = 0, hjust = 0, vjust = 0),
      axis.line.x = element_line(size = 0.3),  
      axis.line.y = element_line(size = 0.3),
      axis.ticks.x = element_line(size = 0.3),
      axis.ticks.y = element_line(size = 0.3)
    ) +
    scale_color_manual(values = colors) +  
    scale_fill_manual(values = colors) + 
    geom_hline(yintercept = 0, linetype = "dashed", size = 0.2, alpha = 1) +
    geom_point(position = position_jitter(0.1), size = 3, alpha = 0.5, aes(shape = 16)) +
    geom_flat_violin(
      position = position_nudge(x = 0.2, y = 0), 
      adjust = 2, 
      alpha = 0.5, 
      trim = FALSE, 
      scale = "width"
    ) +
    geom_boxplot(
      aes(x = as.numeric(Group) + 0.2, y = value), 
      notch = FALSE, 
      lwd = 0.2, 
      width = 0.1, 
      varwidth = FALSE, 
      outlier.shape = NA, 
      alpha = 0.5, 
      colour = "black", 
      show.legend = FALSE
    )
}

#' Process BCR/TCR repertoire files
#'
#' @param file_path Path to directory containing CSV files
#' @param output_file Output filename for combined data
#' @return Combined data frame
#' @export
process_repertoire_files <- function(file_path, output_file) {
  require(readr)
  require(dplyr)
  
  # Get list of CSV files
  files <- list.files(path = file_path, pattern = "*.csv", full.names = TRUE)
  
  # Read and combine files
  combined_df <- bind_rows(lapply(files, function(file) {
    read_csv(file, show_col_types = FALSE)
  }))
  
  # Save combined data
  write_csv(combined_df, output_file)
  
  return(combined_df)
}

#' Calculate isotype percentages for BCR data
#'
#' @param data Data frame with BCR data
#' @param gene_col Column name for gene information
#' @param filter_chains Vector of chain prefixes to filter out
#' @return Data frame with percentages
#' @export
calculate_isotype_percentages <- function(data, gene_col = "c_gene", 
                                          filter_chains = c("IGL", "IGK")) {
  require(dplyr)
  
  # Calculate percentages
  df_percentages <- data %>%
    count(.data[[gene_col]]) %>%
    mutate(percentage = n / sum(n) * 100)
  
  # Filter out specified chains
  for (chain in filter_chains) {
    df_percentages <- df_percentages %>%
      filter(!startsWith(.data[[gene_col]], chain))
  }
  
  # Recalculate percentages after filtering
  df_percentages <- df_percentages %>%
    mutate(
      new_total = sum(n),
      new_percentage = (n / new_total) * 100
    ) %>%
    select(-new_total)
  
  df_percentages$new_percentage <- round(df_percentages$new_percentage, 2)
  
  return(df_percentages)
}

#' Create donut chart for isotype distribution
#'
#' @param data Data frame with isotype percentages
#' @param gene_col Column name for genes/isotypes
#' @param percentage_col Column name for percentages
#' @param title Plot title
#' @param show_labels Logical, whether to show percentage labels
#' @return ggplot object
#' @export
create_isotype_donut <- function(data, gene_col, percentage_col, 
                                 title = "Distribution of Isotypes",
                                 show_labels = FALSE) {
  require(ggplot2)
  require(RColorBrewer)
  
  # Define color palette
  set3_colors <- brewer.pal(12, "Set3")
  isotype_colors <- c(
    "IGHA1" = set3_colors[1],
    "IGHA2" = set3_colors[2],
    "IGHD" = set3_colors[3],
    "IGHG1" = set3_colors[4],
    "IGHG2" = set3_colors[5],
    "IGHM" = set3_colors[6],
    "IGHG3" = set3_colors[7],
    "IGHG4" = set3_colors[8]
  )
  
  # Create base plot
  p <- ggplot(data, aes(x = 2, y = .data[[percentage_col]], fill = .data[[gene_col]])) +
    geom_bar(stat = "identity", width = 1, color = "black") +
    coord_polar("y", start = 0) +
    theme_void() +
    labs(title = title, fill = "Isotype") +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5)
    ) +
    xlim(1, 2.5) +
    scale_fill_manual(values = isotype_colors)
  
  # Add labels if requested
  if (show_labels) {
    p <- p + geom_text(
      aes(label = paste0(round(.data[[percentage_col]], 1), "%")), 
      position = position_stack(vjust = 0.5)
    )
  }
  
  return(p)
}

#' Create enhanced volcano plot for pseudobulk analysis
#'
#' @param data Data frame with DEG results
#' @param gene_col Column name for gene symbols
#' @param fc_col Column name for log2 fold change
#' @param pval_col Column name for adjusted p-values
#' @param title Plot title
#' @param fc_cutoff Fold change cutoff (log2 scale)
#' @param p_cutoff P-value cutoff
#' @param n_top Number of top genes to label
#' @return EnhancedVolcano plot
#' @export
create_enhanced_volcano <- function(data, gene_col = "genes", 
                                    fc_col = "log2FoldChange", 
                                    pval_col = "padj",
                                    title = "", 
                                    fc_cutoff = 0.585,
                                    p_cutoff = 0.05,
                                    n_top = 10) {
  require(EnhancedVolcano)
  require(tidyverse)
  
  # Ensure data has required columns
  if (!gene_col %in% colnames(data)) {
    data <- rownames_to_column(data, var = gene_col)
  }
  
  # Define colors
  keyvals <- ifelse(
    data[[fc_col]] <= -fc_cutoff & data[[pval_col]] <= p_cutoff, 
    "dodgerblue3",
    ifelse(
      data[[fc_col]] >= fc_cutoff & data[[pval_col]] <= p_cutoff, 
      "chartreuse3",
      "grey"
    )
  )
  keyvals[is.na(keyvals)] <- "grey"
  names(keyvals)[keyvals == "chartreuse3"] <- "Upregulated"
  names(keyvals)[keyvals == "grey"] <- "Not significant"
  names(keyvals)[keyvals == "dodgerblue3"] <- "Downregulated"
  
  # Get top genes
  top_up <- data %>% 
    filter(.data[[fc_col]] >= fc_cutoff & .data[[pval_col]] <= p_cutoff) %>% 
    arrange(desc(.data[[fc_col]])) %>% 
    head(n_top)
  
  top_down <- data %>% 
    filter(.data[[fc_col]] <= -fc_cutoff & .data[[pval_col]] <= p_cutoff) %>% 
    arrange(.data[[fc_col]]) %>% 
    head(n_top)
  
  top_genes <- c(top_up[[gene_col]], top_down[[gene_col]])
  
  # Create volcano plot
  vp <- EnhancedVolcano(
    data,
    lab = data[[gene_col]],
    x = fc_col,
    y = pval_col,
    pCutoff = p_cutoff,
    FCcutoff = fc_cutoff,
    pointSize = 3,
    xlab = bquote(~Log[2]~ "fold change"),
    selectLab = top_genes,
    boxedLabels = TRUE,
    labSize = 3,
    drawConnectors = TRUE,
    widthConnectors = 0.5,
    lengthConnectors = unit(0.01, "npc"),
    maxoverlapsConnectors = NULL,
    max.overlaps = 50,
    arrowheads = FALSE,
    shape = 16,
    title = title,
    legendPosition = "right",
    legendLabSize = 14,
    colAlpha = 0.6,
    colCustom = keyvals,
    labFace = "bold"
  )
  
  return(vp)
}

#' Perform statistical testing on cell type proportions
#'
#' @param data_wide Wide-format data frame with proportions
#' @param group_col Column name for group variable
#' @return Data frame with fold changes and p-values
#' @export
calculate_proportion_stats <- function(data_wide, group_col = "Group") {
  require(dplyr)
  
  # Split by group
  groups <- unique(data_wide[[group_col]])
  if (length(groups) != 2) {
    stop("This function requires exactly two groups")
  }
  
  group1 <- data_wide[data_wide[[group_col]] == groups[1], ] %>% select(-all_of(group_col))
  group2 <- data_wide[data_wide[[group_col]] == groups[2], ] %>% select(-all_of(group_col))
  
  # Calculate statistics
  results <- data.frame(
    Cluster = colnames(group1), 
    Log2FC = NA, 
    Log10pvalue = NA
  )
  
  for (i in 1:ncol(group1)) {
    x <- group1[, i]
    y <- group2[, i]
    
    # Calculate fold change
    Log2FC <- log2(median(x, na.rm = TRUE) / median(y, na.rm = TRUE))
    
    # Calculate p-value
    pvalue <- tryCatch({
      wilcox.test(x, y)$p.value
    }, error = function(e) {
      NA
    })
    
    Log10pvalue <- -log10(pvalue)
    
    results[i, 2] <- Log2FC
    results[i, 3] <- Log10pvalue
  }
  
  # Remove NA and infinite values
  results <- results[
    complete.cases(results$Log2FC, results$Log10pvalue) & 
      !is.infinite(results$Log2FC) & 
      !is.infinite(results$Log10pvalue), 
  ]
  
  return(results)
}

#' Load and prepare immunarch data
#'
#' @param file_path Path to repertoire data directory
#' @return Immunarch data object
#' @export
load_repertoire_data <- function(file_path) {
  require(immunarch)
  
  immdata <- repLoad(file_path)
  
  return(immdata)
}

#' Create clonotype abundance plot
#'
#' @param immdata Immunarch data object
#' @param top_clones Vector of top clone numbers to display
#' @return Data frame with formatted clonotype proportions
#' @export
format_clonotype_abundance <- function(immdata, top_clones = c(10, 100, 10000)) {
  require(immunarch)
  require(dplyr)
  
  # Calculate clonality
  imm_top <- repClonality(immdata$data, .method = "top", .head = top_clones)
  
  # Format for plotting
  df <- as.data.frame(imm_top)
  colnames(df) <- c("Top10", "Top100", "Others")
  
  # Adjust categories
  df$Top100 <- df$Top100 - df$Top10
  df$Others <- 1 - (df$Top10 + df$Top100)
  
  # Convert to percentages
  df[, c("Top10", "Top100", "Others")] <- df[, c("Top10", "Top100", "Others")] * 100
  
  return(df)
}

################################################################################
# End of functions
################################################################################

