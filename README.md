# immune-cell-profiling-anti-gad65
Immune cell profiling reveals expanded stem cell-like memory T cells in anti-GAD65-associated neurological syndromes

This repository contains all analysis code associated with the preprint **"Immune cell profiling reveals expanded stem cell-like memory T cells in anti-GAD65-associated neurological syndromes"** [web:16][web:20]

## Overview

The code in this repository reproduces the main computational analyses in the manuscript, including: [web:20]
- Preprocessing and quality control of CSF and blood single-cell RNA-seq data. [web:20]
- Clustering and annotation of immune cell populations, with a focus on CD4+ stem cell-like memory T cells (TSCM). [web:20]
- Analysis of T cell clonal expansion and proinflammatory gene expression. [web:20]
- Integration of B cell receptor (BCR) repertoire data and identification of GAD65-reactive B cell clones. [web:20]
- Generation of figures and summary statistics reported in the preprint. [web:20]

## Repository structure

A suggested directory layout is as follows:

```text
.
├── data/                    # Input data pointers or small example data (no raw patient-level data)
├── envs/                    # Environment / environment.yml files
├── notebooks/               # Jupyter or RMarkdown notebooks for exploratory analyses
├── scripts/                 # Standalone scripts for batch processing
│   ├── 01_qc_preprocessing/
│   ├── 02_integration_clustering/
│   ├── 03_t_cell_clonality/
│   ├── 04_bcr_analysis/
│   └── 05_figures_panels/
├── config/                  # Config files (paths, parameters)
├── results/                 # Output (figures, tables, intermediate objects) – typically not versioned
└── README.md
