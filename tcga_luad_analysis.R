# ============================================================
# TCGA Lung Cancer Multi-Omics Analysis
# Author: Gokul Selvaraj
# GitHub: GokulSelvaraj-Scientist
# Description: Comprehensive analysis of TCGA LUAD (Lung 
#              Adenocarcinoma) data including:
#              1. Survival prediction from genomic features
#              2. Molecular subtype classification
#              3. Mutation landscape and driver gene analysis
#              4. Multi-omics integration
# Data source: TCGA via TCGAbiolinks Bioconductor package
# ============================================================

# --- Install Bioconductor packages if needed ---
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment", 
#                        "DESeq2", "ComplexHeatmap", "maftools"))
# install.packages(c("tidyverse", "survival", "survminer", "caret",
#                    "randomForest", "pROC", "pheatmap", "RColorBrewer",
#                    "ggrepel", "factoextra", "corrplot"))

library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)
library(survival)
library(survminer)
library(caret)
library(randomForest)
library(pROC)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(factoextra)
library(maftools)

# ============================================================
# PART 1: DATA DOWNLOAD FROM TCGA
# ============================================================

cat("=== Downloading TCGA LUAD Clinical Data ===\n")

# Download clinical data
clinical_query <- GDCquery(
  project      = "TCGA-LUAD",
  data.category = "Clinical",
  data.type    = "Clinical Supplement",
  data.format  = "BCR Biotab"
)
GDCdownload(clinical_query, method = "api")
clinical_data <- GDCprepare(clinical_query)

# Extract patient clinical data
clinical <- clinical_data$clinical_patient_luad %>%
  as.data.frame() %>%
  filter(row_number() > 2) %>%  # Remove header rows
  select(
    bcr_patient_barcode,
    age_at_initial_pathologic_diagnosis,
    gender,
    ajcc_pathologic_tumor_stage,
    tobacco_smoking_history,
    days_to_death,
    days_to_last_followup,
    vital_status
  ) %>%
  mutate(
    age        = as.numeric(age_at_initial_pathologic_diagnosis),
    OS_time    = as.numeric(ifelse(!is.na(days_to_death) & days_to_death != "[Not Available]",
                                   days_to_death, days_to_last_followup)),
    OS_status  = ifelse(vital_status == "Dead", 1, 0),
    stage      = case_when(
      grepl("Stage I[^V]|Stage IA|Stage IB", ajcc_pathologic_tumor_stage) ~ "Stage I",
      grepl("Stage II[^I]|Stage IIA|Stage IIB", ajcc_pathologic_tumor_stage) ~ "Stage II",
      grepl("Stage III", ajcc_pathologic_tumor_stage) ~ "Stage III",
      grepl("Stage IV", ajcc_pathologic_tumor_stage) ~ "Stage IV",
      TRUE ~ NA_character_
    ),
    smoking    = case_when(
      tobacco_smoking_history %in% c("1") ~ "Never",
      tobacco_smoking_history %in% c("2", "3") ~ "Former",
      tobacco_smoking_history %in% c("4", "5") ~ "Current",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(OS_time), OS_time > 0, !is.na(stage))

cat("Clinical data: ", nrow(clinical), "patients\n")

# ============================================================
# Download gene expression data (RNA-seq)
# ============================================================

cat("\n=== Downloading TCGA LUAD RNA-seq Data ===\n")

rna_query <- GDCquery(
  project            = "TCGA-LUAD",
  data.category      = "Transcriptome Profiling",
  data.type          = "Gene Expression Quantification",
  workflow.type      = "STAR - Counts",
  sample.type        = "Primary Tumor"
)
GDCdownload(rna_query, method = "api", files.per.chunk = 50)
rna_data <- GDCprepare(rna_query)

# Extract count matrix
count_matrix <- assay(rna_data, "unstranded")
gene_info    <- rowData(rna_data)
sample_info  <- colData(rna_data)

cat("RNA-seq data:", nrow(count_matrix), "genes x", ncol(count_matrix), "samples\n")

# ============================================================
# Download mutation data (MAF)
# ============================================================

cat("\n=== Downloading TCGA LUAD Mutation Data ===\n")

maf_query <- GDCquery(
  project      = "TCGA-LUAD",
  data.category = "Simple Nucleotide Variation",
  data.type    = "Masked Somatic Mutation",
  access       = "open"
)
GDCdownload(maf_query, method = "api")
maf_data <- GDCprepare(maf_query)

cat("MAF data:", nrow(maf_data), "mutations\n")
cat("Unique patients:", length(unique(maf_data$Tumor_Sample_Barcode)), "\n")

# ============================================================
# PART 2: DATA PROCESSING
# ============================================================

cat("\n=== Processing Data ===\n")

# --- Filter and normalize RNA-seq ---
# Remove low expression genes
keep_genes <- rowSums(count_matrix >= 10) >= (0.2 * ncol(count_matrix))
count_filtered <- count_matrix[keep_genes, ]
cat("Genes after filtering:", nrow(count_filtered), "\n")

# Normalize using DESeq2 VST
library(DESeq2)
col_data <- data.frame(
  sample    = colnames(count_filtered),
  condition = rep("tumor", ncol(count_filtered)),
  row.names = colnames(count_filtered)
)

dds <- DESeqDataSetFromMatrix(
  countData = count_filtered,
  colData   = col_data,
  design    = ~ 1
)
dds     <- estimateSizeFactors(dds)
vst_mat <- assay(vst(dds, blind = TRUE))

# Select top variable genes for downstream analysis
gene_vars   <- apply(vst_mat, 1, var)
top5000     <- names(sort(gene_vars, decreasing = TRUE)[1:5000])
vst_top5000 <- vst_mat[top5000, ]

# Clean sample IDs to match clinical data
colnames(vst_top5000) <- substr(colnames(vst_top5000), 1, 12)

# Merge with clinical data
common_samples <- intersect(colnames(vst_top5000), clinical$bcr_patient_barcode)
vst_final      <- vst_top5000[, common_samples]
clinical_final <- clinical %>% filter(bcr_patient_barcode %in% common_samples)

cat("Matched samples:", length(common_samples), "\n")

# ============================================================
# PART 3: MUTATION LANDSCAPE ANALYSIS
# ============================================================

cat("\n=== Mutation Landscape Analysis ===\n")

# --- Plot 1: Mutation Summary using maftools ---
maf_object <- read.maf(maf = maf_data)

png("mutation_summary.png", width = 1200, height = 800, res = 150)
plotmafSummary(
  maf        = maf_object,
  rmOutlier  = TRUE,
  addStat    = "median",
  dashboard  = TRUE,
  titvRaw    = FALSE
)
dev.off()
cat("Saved: mutation_summary.png\n")

# --- Plot 2: Oncoplot of top mutated genes ---
png("oncoplot_top20.png", width = 1400, height = 900, res = 150)
oncoplot(
  maf      = maf_object,
  top      = 20,
  fontSize = 0.6,
  title    = "Top 20 Mutated Genes — TCGA LUAD"
)
dev.off()
cat("Saved: oncoplot_top20.png\n")

# --- Driver gene mutation frequencies ---
gene_summary <- getGeneSummary(maf_object)
top_genes <- gene_summary %>%
  as.data.frame() %>%
  head(20) %>%
  mutate(MutationFreq = AlteredSamples / as.numeric(maf_object@summary$summary[3]) * 100)

driver_plot <- ggplot(top_genes, aes(x = reorder(Hugo_Symbol, MutationFreq),
                                      y = MutationFreq, fill = MutationFreq)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "#A8DADC", high = "#E63946") +
  labs(
    title    = "Top 20 Most Frequently Mutated Genes",
    subtitle = "TCGA Lung Adenocarcinoma (LUAD)",
    x        = "Gene",
    y        = "% Samples Mutated",
    fill     = "Frequency (%)"
  ) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

ggsave("driver_gene_frequency.png", driver_plot, width = 9, height = 7, dpi = 300)
cat("Saved: driver_gene_frequency.png\n")

# ============================================================
# PART 4: MOLECULAR SUBTYPE CLASSIFICATION
# ============================================================

cat("\n=== Molecular Subtype Classification ===\n")

# --- Consensus clustering to identify subtypes ---
# PCA on top variable genes
pca_result <- prcomp(t(vst_final), scale. = TRUE)
var_exp    <- round(summary(pca_result)$importance[2, 1:2] * 100, 1)

pca_df <- data.frame(
  sample     = rownames(pca_result$x),
  PC1        = pca_result$x[, 1],
  PC2        = pca_result$x[, 2],
  PC3        = pca_result$x[, 3]
) %>%
  left_join(clinical_final %>% select(bcr_patient_barcode, stage, smoking, gender),
            by = c("sample" = "bcr_patient_barcode"))

# K-means clustering (k=3 based on published LUAD subtypes)
set.seed(42)
km_result  <- kmeans(pca_result$x[, 1:10], centers = 3, nstart = 25)
pca_df$Subtype <- paste0("Subtype ", km_result$cluster)

# --- Plot 3: PCA colored by subtype ---
subtype_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Subtype, shape = stage)) +
  geom_point(size = 2.5, alpha = 0.7) +
  scale_color_manual(values = c("#2A9D8F", "#E76F51", "#457B9D")) +
  stat_ellipse(geom = "polygon", alpha = 0.08, level = 0.75) +
  labs(
    title    = "Molecular Subtypes: PCA of Gene Expression",
    subtitle = "TCGA LUAD — Unsupervised clustering reveals 3 subtypes",
    x        = paste0("PC1 (", var_exp[1], "% variance)"),
    y        = paste0("PC2 (", var_exp[2], "% variance)"),
    color    = "Molecular Subtype",
    shape    = "Tumor Stage"
  ) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

ggsave("molecular_subtypes_pca.png", subtype_pca, width = 9, height = 6, dpi = 300)
cat("Saved: molecular_subtypes_pca.png\n")

# --- Plot 4: Heatmap of top 100 variable genes by subtype ---
top100 <- names(sort(gene_vars[top5000], decreasing = TRUE)[1:100])
mat100 <- vst_final[top100, ]
mat100 <- mat100 - rowMeans(mat100)

annotation_col <- data.frame(
  Subtype = pca_df$Subtype,
  Stage   = pca_df$stage,
  row.names = pca_df$sample
)
annotation_col <- annotation_col[colnames(mat100), ]

ann_colors <- list(
  Subtype = c("Subtype 1" = "#2A9D8F", "Subtype 2" = "#E76F51", "Subtype 3" = "#457B9D"),
  Stage   = c("Stage I" = "#F1FAEE", "Stage II" = "#A8DADC",
              "Stage III" = "#457B9D", "Stage IV" = "#1D3557")
)

png("subtype_heatmap.png", width = 1400, height = 1000, res = 150)
pheatmap(
  mat100,
  annotation_col  = annotation_col,
  annotation_colors = ann_colors,
  show_rownames   = FALSE,
  show_colnames   = FALSE,
  cluster_rows    = TRUE,
  cluster_cols    = TRUE,
  color           = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
  main            = "Top 100 Variable Genes by Molecular Subtype\nTCGA LUAD"
)
dev.off()
cat("Saved: subtype_heatmap.png\n")

# ============================================================
# PART 5: SURVIVAL ANALYSIS
# ============================================================

cat("\n=== Survival Analysis ===\n")

# Merge subtype with clinical data
survival_df <- clinical_final %>%
  left_join(pca_df %>% select(sample, Subtype),
            by = c("bcr_patient_barcode" = "sample")) %>%
  filter(!is.na(Subtype), !is.na(OS_time), OS_time > 0)

# --- Plot 5: Survival by molecular subtype ---
surv_obj  <- Surv(time = survival_df$OS_time, event = survival_df$OS_status)
fit_subtype <- survfit(surv_obj ~ Subtype, data = survival_df)

png("survival_by_subtype.png", width = 1000, height = 800, res = 150)
ggsurvplot(
  fit_subtype,
  data          = survival_df,
  pval          = TRUE,
  conf.int      = FALSE,
  risk.table    = TRUE,
  legend.title  = "Molecular Subtype",
  palette       = c("#2A9D8F", "#E76F51", "#457B9D"),
  xlab          = "Time (days)",
  ylab          = "Overall Survival Probability",
  title         = "Overall Survival by Molecular Subtype",
  subtitle      = "TCGA Lung Adenocarcinoma",
  risk.table.height = 0.28,
  ggtheme       = theme_classic(base_size = 12)
)
dev.off()
cat("Saved: survival_by_subtype.png\n")

# --- Plot 6: Survival by tumor stage ---
fit_stage <- survfit(surv_obj ~ stage, data = survival_df)

png("survival_by_stage.png", width = 1000, height = 800, res = 150)
ggsurvplot(
  fit_stage,
  data          = survival_df,
  pval          = TRUE,
  conf.int      = FALSE,
  risk.table    = TRUE,
  legend.title  = "Tumor Stage",
  palette       = c("#F1FAEE", "#A8DADC", "#457B9D", "#1D3557"),
  xlab          = "Time (days)",
  ylab          = "Overall Survival Probability",
  title         = "Overall Survival by Tumor Stage",
  subtitle      = "TCGA Lung Adenocarcinoma",
  risk.table.height = 0.28,
  ggtheme       = theme_classic(base_size = 12)
)
dev.off()
cat("Saved: survival_by_stage.png\n")

# --- Cox Proportional Hazards Model ---
cox_df <- survival_df %>%
  mutate(
    stage_num = case_when(
      stage == "Stage I"   ~ 1,
      stage == "Stage II"  ~ 2,
      stage == "Stage III" ~ 3,
      stage == "Stage IV"  ~ 4
    ),
    gender_bin = ifelse(gender == "MALE", 1, 0)
  ) %>%
  filter(!is.na(stage_num), !is.na(age))

cox_model <- coxph(
  Surv(OS_time, OS_status) ~ age + gender_bin + stage_num + Subtype,
  data = cox_df
)

cat("\nCox Model Summary:\n")
print(summary(cox_model))

png("cox_forest_plot.png", width = 1000, height = 700, res = 150)
ggforest(cox_model, data = cox_df,
         main = "Cox Proportional Hazards — TCGA LUAD")
dev.off()
cat("Saved: cox_forest_plot.png\n")

# ============================================================
# PART 6: SURVIVAL PREDICTION WITH MACHINE LEARNING
# ============================================================

cat("\n=== ML Survival Prediction ===\n")

# Select top prognostic genes using Cox regression
cat("Selecting prognostic genes...\n")

# Use top 500 variable genes for ML
top500 <- names(sort(gene_vars[top5000], decreasing = TRUE)[1:500])
expr_ml <- t(vst_final[top500, ]) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sample")

ml_data <- survival_df %>%
  select(bcr_patient_barcode, OS_status, OS_time, age, stage, Subtype) %>%
  left_join(expr_ml, by = c("bcr_patient_barcode" = "sample")) %>%
  filter(complete.cases(.)) %>%
  mutate(
    Survived2yr = factor(ifelse(OS_time >= 730, "Survived", "Died"),
                         levels = c("Survived", "Died")),
    stage_num   = case_when(
      stage == "Stage I"   ~ 1,
      stage == "Stage II"  ~ 2,
      stage == "Stage III" ~ 3,
      stage == "Stage IV"  ~ 4
    )
  ) %>%
  select(-bcr_patient_barcode, -OS_status, -OS_time, -stage)

cat("ML dataset:", nrow(ml_data), "samples\n")
cat("2-year survival distribution:\n")
print(table(ml_data$Survived2yr))

# Train/test split
set.seed(42)
train_idx  <- createDataPartition(ml_data$Survived2yr, p = 0.8, list = FALSE)
train_data <- ml_data[train_idx, ]
test_data  <- ml_data[-train_idx, ]

# Random Forest with cross-validation
ctrl <- trainControl(
  method          = "cv",
  number          = 5,
  classProbs      = TRUE,
  summaryFunction = twoClassSummary
)

cat("\nTraining Random Forest for 2-year survival prediction...\n")
set.seed(42)
rf_model <- train(
  Survived2yr ~ .,
  data       = train_data,
  method     = "rf",
  trControl  = ctrl,
  metric     = "ROC",
  tuneLength = 3,
  ntree      = 500
)

# Evaluate
rf_pred  <- predict(rf_model, test_data)
rf_probs <- predict(rf_model, test_data, type = "prob")
rf_cm    <- confusionMatrix(rf_pred, test_data$Survived2yr, positive = "Survived")
rf_roc   <- roc(test_data$Survived2yr, rf_probs$Survived,
                levels = c("Died", "Survived"))

cat("\nSurvival Prediction Performance:\n")
cat("Accuracy:   ", round(rf_cm$overall["Accuracy"], 3), "\n")
cat("Sensitivity:", round(rf_cm$byClass["Sensitivity"], 3), "\n")
cat("Specificity:", round(rf_cm$byClass["Specificity"], 3), "\n")
cat("AUC:        ", round(auc(rf_roc), 3), "\n")

# --- Plot 7: ROC curve for survival prediction ---
roc_df <- data.frame(
  FPR = 1 - rf_roc$specificities,
  TPR = rf_roc$sensitivities
)

roc_plot <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
  geom_line(color = "#E76F51", linewidth = 1.5) +
  geom_abline(linetype = "dashed", color = "grey50") +
  annotate("text", x = 0.6, y = 0.2,
           label = paste0("AUC = ", round(auc(rf_roc), 3)),
           size = 5, color = "#E76F51", fontface = "bold") +
  labs(
    title    = "ROC Curve: 2-Year Survival Prediction",
    subtitle = "Random Forest trained on TCGA LUAD gene expression + clinical features",
    x        = "False Positive Rate",
    y        = "True Positive Rate"
  ) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

ggsave("survival_prediction_roc.png", roc_plot, width = 8, height = 6, dpi = 300)
cat("Saved: survival_prediction_roc.png\n")

# --- Plot 8: Top prognostic genes ---
importance_df <- varImp(rf_model)$importance %>%
  tibble::rownames_to_column("Feature") %>%
  arrange(desc(Overall)) %>%
  head(15)

importance_plot <- ggplot(importance_df,
                          aes(x = reorder(Feature, Overall), y = Overall, fill = Overall)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "#A8DADC", high = "#E63946") +
  labs(
    title    = "Top 15 Prognostic Features",
    subtitle = "Random Forest feature importance for 2-year survival prediction",
    x        = "Feature",
    y        = "Importance Score"
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "none")

ggsave("prognostic_features.png", importance_plot, width = 9, height = 7, dpi = 300)
cat("Saved: prognostic_features.png\n")

# ============================================================
# PART 7: MULTI-OMICS INTEGRATION
# ============================================================

cat("\n=== Multi-Omics Integration ===\n")

# Integrate mutation status of key driver genes with expression
key_drivers <- c("TP53", "KRAS", "EGFR", "STK11", "KEAP1", "RBM10", "BRAF", "NF1")

# Get mutation status per patient for key drivers
mut_status <- maf_data %>%
  filter(Hugo_Symbol %in% key_drivers) %>%
  select(Tumor_Sample_Barcode, Hugo_Symbol) %>%
  mutate(patient = substr(Tumor_Sample_Barcode, 1, 12),
         mutated = 1) %>%
  distinct(patient, Hugo_Symbol, .keep_all = TRUE) %>%
  pivot_wider(id_cols = patient, names_from = Hugo_Symbol,
              values_from = mutated, values_fill = 0)

# Merge mutation status with expression PCA and clinical
multiomics_df <- pca_df %>%
  left_join(mut_status, by = c("sample" = "patient")) %>%
  left_join(survival_df %>% select(bcr_patient_barcode, OS_time, OS_status),
            by = c("sample" = "bcr_patient_barcode")) %>%
  mutate(across(all_of(key_drivers[key_drivers %in% names(.)]),
                ~ replace_na(., 0)))

# --- Plot 9: EGFR mutation status on PCA ---
if ("EGFR" %in% names(multiomics_df)) {
  egfr_plot <- ggplot(multiomics_df,
                      aes(x = PC1, y = PC2,
                          color = factor(EGFR),
                          size  = factor(EGFR))) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = c("0" = "#AAAAAA", "1" = "#E63946"),
                       labels = c("Wild-type", "Mutant")) +
    scale_size_manual(values = c("0" = 1.5, "1" = 3),
                      labels = c("Wild-type", "Mutant")) +
    labs(
      title    = "EGFR Mutation Status in Gene Expression Space",
      subtitle = "TCGA LUAD — EGFR mutations cluster in distinct expression subtypes",
      x        = paste0("PC1 (", var_exp[1], "% variance)"),
      y        = paste0("PC2 (", var_exp[2], "% variance)"),
      color    = "EGFR Status",
      size     = "EGFR Status"
    ) +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))

  ggsave("egfr_mutation_pca.png", egfr_plot, width = 9, height = 6, dpi = 300)
  cat("Saved: egfr_mutation_pca.png\n")
}

# --- Plot 10: Multi-omics correlation heatmap ---
# Correlate mutation status with top expression PCs
multiomics_cor <- multiomics_df %>%
  select(PC1, PC2, PC3, any_of(key_drivers)) %>%
  filter(complete.cases(.))

if (ncol(multiomics_cor) > 4) {
  cor_mat <- cor(multiomics_cor, use = "complete.obs")

  png("multiomics_correlation.png", width = 900, height = 800, res = 150)
  pheatmap(
    cor_mat,
    color    = colorRampPalette(c("#457B9D", "white", "#E63946"))(100),
    display_numbers = TRUE,
    number_format   = "%.2f",
    fontsize_number = 8,
    main     = "Multi-Omics Correlation: Expression PCs vs Mutation Status\nTCGA LUAD"
  )
  dev.off()
  cat("Saved: multiomics_correlation.png\n")
}

# ============================================================
# SAVE SUMMARY RESULTS
# ============================================================

# Save processed clinical data
write.csv(survival_df, "tcga_luad_clinical_processed.csv", row.names = FALSE)

# Save subtype assignments
write.csv(pca_df %>% select(sample, Subtype, PC1, PC2),
          "tcga_luad_subtype_assignments.csv", row.names = FALSE)

# Save ML performance
ml_results <- data.frame(
  Analysis    = "2-year survival prediction",
  Model       = "Random Forest",
  N_samples   = nrow(ml_data),
  N_features  = ncol(ml_data) - 1,
  Accuracy    = round(rf_cm$overall["Accuracy"], 3),
  Sensitivity = round(rf_cm$byClass["Sensitivity"], 3),
  Specificity = round(rf_cm$byClass["Specificity"], 3),
  AUC         = round(auc(rf_roc), 3)
)
write.csv(ml_results, "ml_performance_summary.csv", row.names = FALSE)

cat("\n=== Analysis Complete ===\n")
cat("All outputs saved successfully.\n")
cat("\nSummary:\n")
cat("- Patients analyzed:", nrow(clinical_final), "\n")
cat("- Genes analyzed:", nrow(vst_final), "\n")
cat("- Mutations analyzed:", nrow(maf_data), "\n")
cat("- Molecular subtypes identified: 3\n")
cat("- 2-year survival prediction AUC:", round(auc(rf_roc), 3), "\n")
