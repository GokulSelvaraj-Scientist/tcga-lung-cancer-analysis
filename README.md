# TCGA Lung Adenocarcinoma: Multi-Omics Analysis

## Overview
A comprehensive multi-omics analysis of lung adenocarcinoma (LUAD) using real patient data from The Cancer Genome Atlas (TCGA). This project integrates gene expression, somatic mutation, and clinical data from hundreds of real lung cancer patients to characterize molecular subtypes, identify driver genes, predict patient survival, and explore multi-omics relationships.

## Why This Matters
Lung adenocarcinoma is the most common subtype of non-small cell lung cancer (NSCLC) and the leading cause of cancer mortality worldwide. Despite significant advances with targeted therapies (EGFR inhibitors, ALK inhibitors) and immunotherapy (checkpoint inhibitors), most patients develop resistance and prognosis remains poor. Multi-omics analysis of TCGA data is a cornerstone approach for:

- **Biomarker discovery** — identifying molecular features associated with survival enables development of companion diagnostics for patient selection in clinical trials
- **Molecular subtyping** — discovering transcriptionally distinct patient subgroups helps explain heterogeneous treatment responses and may define new indications for targeted agents
- **Driver gene characterization** — understanding mutation frequencies and co-occurrence patterns across a large real-world cohort informs drug target prioritization and combination therapy rationale
- **Survival prediction** — machine learning models trained on real TCGA data provide proof-of-concept for prognostic signatures that can be prospectively validated
- **Precision oncology** — integrating genomic and transcriptomic data mirrors the analytical framework used in FDA submissions for companion diagnostics and tumor-agnostic approvals

This is the type of analysis performed in translational research groups at major cancer centers and pharmaceutical companies working on oncology drug development programs.

## Data Source
- **Database:** The Cancer Genome Atlas (TCGA) — accessed via TCGAbiolinks Bioconductor package
- **Project:** TCGA-LUAD (Lung Adenocarcinoma)
- **Data types:** RNA-seq gene expression, somatic mutations (MAF), clinical data
- **Access:** Publicly available via GDC Data Portal (https://portal.gdc.cancer.gov/)

## Methods

### 1. Mutation Landscape Analysis
- Somatic mutation profiling using maftools
- Mutation summary (SNV types, TMB distribution)
- Oncoplot of top 20 mutated genes
- Driver gene mutation frequency analysis

### 2. Molecular Subtype Classification
- Variance-based gene filtering (top 5000 variable genes)
- DESeq2 variance stabilizing transformation (VST)
- PCA dimensionality reduction
- K-means consensus clustering (k=3 subtypes)
- Subtype characterization by heatmap

### 3. Survival Analysis
- Kaplan-Meier survival curves by molecular subtype and tumor stage
- Log-rank tests for group comparisons
- Cox proportional hazards regression adjusting for age, sex, stage, and subtype

### 4. Survival Prediction with Machine Learning
- Binary classification: 2-year overall survival (survived vs died)
- Features: top 500 variable genes + clinical variables
- Random Forest with 5-fold cross-validation
- Performance: Accuracy, Sensitivity, Specificity, AUC-ROC

### 5. Multi-Omics Integration
- Overlay of somatic mutation status on gene expression PCA space
- Correlation analysis between expression principal components and mutation status of key driver genes
- Visualization of EGFR mutation distribution across molecular subtypes

## Outputs
| File | Description |
|---|---|
| `mutation_summary.png` | Mutation landscape summary (TMB, variant types) |
| `oncoplot_top20.png` | Oncoplot of top 20 mutated genes |
| `driver_gene_frequency.png` | Mutation frequency barplot for driver genes |
| `molecular_subtypes_pca.png` | PCA plot colored by molecular subtype |
| `subtype_heatmap.png` | Heatmap of top 100 variable genes by subtype |
| `survival_by_subtype.png` | KM curves by molecular subtype |
| `survival_by_stage.png` | KM curves by tumor stage |
| `cox_forest_plot.png` | Cox model hazard ratios |
| `survival_prediction_roc.png` | ROC curve for ML survival prediction |
| `prognostic_features.png` | Top 15 prognostic gene features |
| `egfr_mutation_pca.png` | EGFR mutation status in expression space |
| `multiomics_correlation.png` | Multi-omics correlation heatmap |
| `tcga_luad_clinical_processed.csv` | Processed clinical data |
| `tcga_luad_subtype_assignments.csv` | Molecular subtype assignments |
| `ml_performance_summary.csv` | ML model performance metrics |

## Key Expected Findings
- **TP53, KRAS, EGFR, STK11** are among the most frequently mutated genes in LUAD — consistent with published literature
- **EGFR mutations** cluster in specific molecular subtypes — reflecting the transcriptional consequences of EGFR-driven oncogenesis
- **Molecular subtypes** show distinct survival outcomes — supporting their clinical relevance as prognostic biomarkers
- **Tumor stage** is a strong independent predictor of survival — Stage I patients have significantly better outcomes than Stage III/IV
- **Gene expression features** improve survival prediction beyond clinical variables alone — supporting the development of transcriptomic prognostic signatures

## Biological and Clinical Interpretation
The integration of mutation and expression data reveals how somatic alterations in key driver genes shape the transcriptional landscape of lung tumors. EGFR mutations — the target of first-line therapies like erlotinib and osimertinib — cluster in distinct expression subtypes, suggesting that transcriptomic profiling can serve as a surrogate for mutation status in some contexts. The identification of 3 molecular subtypes with distinct survival outcomes provides a biological basis for patient stratification in clinical trials, and the prognostic gene signature derived from TCGA data represents a candidate biomarker panel for prospective validation.

## How to Run

### Step 1 — Install packages
```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment", "DESeq2", 
                       "ComplexHeatmap", "maftools"))
install.packages(c("tidyverse", "survival", "survminer", "caret",
                   "randomForest", "pROC", "pheatmap", "RColorBrewer",
                   "ggrepel", "factoextra", "corrplot"))
```

### Step 2 — Run analysis
```r
setwd("path/to/tcga-lung-cancer-analysis")
source("tcga_luad_analysis.R")
```

**Note:** Data download requires internet connection and approximately 2-4 GB of disk space. Initial download takes 30-60 minutes depending on connection speed. Downloaded data is cached locally for subsequent runs.

## Requirements
- R >= 4.0
- TCGAbiolinks, SummarizedExperiment, DESeq2, maftools
- tidyverse, survival, survminer, caret, randomForest, pROC
- pheatmap, RColorBrewer, ggrepel, factoextra

## Author
**Gokul Selvaraj, PhD**
GitHub: [GokulSelvaraj-Scientist](https://github.com/GokulSelvaraj-Scientist)
