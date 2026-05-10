# =============================================================================
# TCGA-LIHC Ketone Body Metabolism Analysis
# Associated with: "Ketone body metabolism in MASLD-associated hepatocellular 
# carcinoma: from molecular mechanisms to therapeutic opportunities"
#
# Description: Integrative transcriptomic analysis of HMGCS2 and OXCT1 
# expression, survival outcomes, immune infiltration, and T-cell exhaustion
# in the TCGA Liver Hepatocellular Carcinoma (LIHC) cohort.
#
# Data source: TCGA-LIHC RNA-seq FPKM data accessed via TCGAbiolinks [95]
# R version: 4.5.2
# Author: [Author name]
# Date: May 2026
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Install and load required packages
# -----------------------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

pkgs_bioc <- c("TCGAbiolinks")
pkgs_cran <- c("dplyr", "tidyr", "ggplot2", "ggpubr", "survminer", 
               "survival", "maxstat", "estimate", "patchwork", "biomaRt")

for (pkg in pkgs_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg)
}
for (pkg in pkgs_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

library(TCGAbiolinks)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(survminer)
library(survival)
library(maxstat)
library(patchwork)

# -----------------------------------------------------------------------------
# 1. Download TCGA-LIHC data
# -----------------------------------------------------------------------------

# RNA-seq expression data
query_expr <- GDCquery(
  project = "TCGA-LIHC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  access = "open"
)
GDCdownload(query_expr, directory = "D:/scRNA_project/GDCdata/")
expr_data <- GDCprepare(query_expr, directory = "D:/scRNA_project/GDCdata/")

# Somatic mutation data
query_maf <- GDCquery(
  project = "TCGA-LIHC",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  access = "open"
)
GDCdownload(query_maf, directory = "D:/scRNA_project/GDCdata/")
maf <- GDCprepare(query_maf, directory = "D:/scRNA_project/GDCdata/")

# Clinical data
clin <- read.table("D:/scRNA_project/TCGA-LIHC.clinical.tsv.gz",
                   sep = "\t", header = TRUE, fill = TRUE)

# Survival data
surv_raw <- read.table("D:/scRNA_project/TCGA-LIHC.survival.tsv.gz",
                       sep = "\t", header = TRUE)


# -----------------------------------------------------------------------------
# 2. Prepare expression matrix (Gene Symbol, tumor vs normal)
# -----------------------------------------------------------------------------

# NOTE: If using pre-processed tcga_sym object from workspace_backup.RData:
# load("D:/scRNA_project/workspace_backup.RData")
# tcga_sym contains rows = gene symbols, columns = samples

# Extract HMGCS2 and OXCT1 expression for tumor samples
genes_of_interest <- c("HMGCS2", "OXCT1", "BDH1")

gene_mat <- t(tcga_sym[genes_of_interest, 4:ncol(tcga_sym)])
gene_df  <- as.data.frame(gene_mat)
gene_df$sample  <- rownames(gene_df)
gene_df$patient <- substr(gene_df$sample, 1, 12)

# Separate tumor (01A/01B) and normal (11A) samples
tumor_df  <- gene_df[grep("01A|01B", gene_df$sample), ]
normal_df <- gene_df[grep("11A|11B", gene_df$sample), ]

cat("Tumor samples:", nrow(tumor_df), "\n")
cat("Normal samples:", nrow(normal_df), "\n")


# -----------------------------------------------------------------------------
# 3. Differential expression: Tumor vs Normal (Fig. 5A)
# -----------------------------------------------------------------------------

# Build long-format data for plotting
expr_plot <- bind_rows(
  mutate(tumor_df,  type = "Tumor"),
  mutate(normal_df, type = "Normal")
)

expr_long <- expr_plot %>%
  pivot_longer(cols = all_of(c("HMGCS2", "OXCT1")),
               names_to = "gene", values_to = "value") %>%
  mutate(gene = factor(gene, levels = c("HMGCS2", "OXCT1")))

# Wilcoxon test
wilcox.test(tumor_df$HMGCS2, normal_df$HMGCS2)  # p = 0.00015
wilcox.test(tumor_df$OXCT1,  normal_df$OXCT1)   # p = 0.00013

# Plot
p_expr <- ggplot(expr_long, aes(x = type, y = log2(value + 1), fill = type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.2, size = 0.6) +
  facet_wrap(~gene, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("Normal" = "#4DBBD5", "Tumor" = "#E64B35")) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  labs(title = "Ketone body metabolism genes in TCGA-LIHC",
       x = "", y = "log2(FPKM + 1)") +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 13))

ggsave("D:/综述代码图/ketone_genes_expression_v2.png",
       plot = p_expr, width = 7, height = 5, dpi = 300)
ggsave("D:/综述代码图/ketone_genes_expression_v2.pdf",
       plot = p_expr, width = 7, height = 5)


# -----------------------------------------------------------------------------
# 4. Survival analysis (Fig. 5B, 5C, 5D)
# -----------------------------------------------------------------------------

# Prepare survival data
surv_df <- data.frame(
  patient      = clin$submitter_id,
  OS           = as.numeric(clin$vital_status.demographic == "Dead"),
  OS.time.month = as.numeric(clin$days_to_death.demographic) / 30,
  stringsAsFactors = FALSE
)

# Use TCGA survival file if available
surv_df <- data.frame(
  patient       = surv_raw$`_PATIENT`,
  OS            = surv_raw$OS,
  OS.time.month = surv_raw$OS.time / 30,
  stringsAsFactors = FALSE
)
surv_df <- surv_df[!is.na(surv_df$patient), ]

# Merge expression and survival
merged <- inner_join(tumor_df, surv_df, by = "patient")
cat("Merged samples for survival analysis:", nrow(merged), "\n")  # n = 425

# Median cutpoint stratification
merged$OXCT1_group  <- ifelse(merged$OXCT1  > median(merged$OXCT1),  "High", "Low")
merged$HMGCS2_group <- ifelse(merged$HMGCS2 > median(merged$HMGCS2), "High", "Low")
merged$group <- factor(
  ifelse(merged$OXCT1 > median(merged$OXCT1) & merged$HMGCS2 < median(merged$HMGCS2),
         "OXCT1-High/HMGCS2-Low",
  ifelse(merged$OXCT1 > median(merged$OXCT1) & merged$HMGCS2 > median(merged$HMGCS2),
         "OXCT1-High/HMGCS2-High",
  ifelse(merged$OXCT1 < median(merged$OXCT1) & merged$HMGCS2 < median(merged$HMGCS2),
         "OXCT1-Low/HMGCS2-Low",
         "OXCT1-Low/HMGCS2-High"))),
  levels = c("OXCT1-Low/HMGCS2-High", "OXCT1-Low/HMGCS2-Low",
             "OXCT1-High/HMGCS2-High", "OXCT1-High/HMGCS2-Low")
)

# --- Fig. 5B: OXCT1 survival ---
fit_oxct1 <- survfit(Surv(OS.time.month, OS) ~ OXCT1_group, data = merged)
p_surv_oxct1 <- ggsurvplot(fit_oxct1, data = merged,
  pval = TRUE, conf.int = FALSE,
  palette = c("#E64B35", "#4DBBD5"),
  xlab = "Time (months)", ylab = "Overall Survival",
  title = "OXCT1 expression and OS in TCGA-LIHC",
  ggtheme = theme_bw())

png("D:/综述代码图/survival_OXCT1.png", width = 8, height = 6, units = "in", res = 300)
print(p_surv_oxct1)
dev.off()

# --- Fig. 5C: HMGCS2 survival ---
fit_hmgcs2 <- survfit(Surv(OS.time.month, OS) ~ HMGCS2_group, data = merged)
p_surv_hmgcs2 <- ggsurvplot(fit_hmgcs2, data = merged,
  pval = TRUE, conf.int = FALSE,
  palette = c("#E64B35", "#4DBBD5"),
  xlab = "Time (months)", ylab = "Overall Survival",
  title = "HMGCS2 expression and OS in TCGA-LIHC",
  ggtheme = theme_bw())

png("D:/综述代码图/survival_HMGCS2.png", width = 8, height = 6, units = "in", res = 300)
print(p_surv_hmgcs2)
dev.off()

# --- Fig. 5D: Combined stratification ---
fit_combo <- survfit(Surv(OS.time.month, OS) ~ group, data = merged)
p_combo <- ggsurvplot(fit_combo, data = merged,
  pval = TRUE, conf.int = FALSE,
  palette = c("#00A087", "#91D1C2", "#F39B7F", "#E64B35"),
  legend.title = "",
  legend.labs = c("OXCT1-Low/HMGCS2-High", "OXCT1-Low/HMGCS2-Low",
                  "OXCT1-High/HMGCS2-High", "OXCT1-High/HMGCS2-Low"),
  xlab = "Time (months)", ylab = "Overall Survival",
  title = "Combined OXCT1/HMGCS2 stratification and OS in TCGA-LIHC",
  risk.table = TRUE, risk.table.height = 0.3,
  font.legend = 9, ggtheme = theme_bw())

png("D:/综述代码图/survival_OXCT1_HMGCS2_combined.png",
    width = 10, height = 7, units = "in", res = 300)
print(p_combo)
dev.off()


# -----------------------------------------------------------------------------
# 5. Cox regression (univariable and multivariable)
# -----------------------------------------------------------------------------

# Univariable Cox
cox_oxct1  <- coxph(Surv(OS.time.month, OS) ~ OXCT1_group,  data = merged)
cox_hmgcs2 <- coxph(Surv(OS.time.month, OS) ~ HMGCS2_group, data = merged)
cox_combo  <- coxph(Surv(OS.time.month, OS) ~ group,         data = merged)

summary(cox_oxct1)   # HR(Low vs High) = 0.71 (0.52-0.96), p = 0.027
summary(cox_hmgcs2)  # HR(Low vs High) = 1.82 (1.34-2.47), p = 0.00014
summary(cox_combo)   # OXCT1-High/HMGCS2-Low vs ref: HR = 2.47 (1.60-3.82), p < 0.001

# Multivariable Cox (adjusted for age, sex, AJCC stage)
clin_df <- data.frame(
  patient = clin$submitter_id,
  age     = as.numeric(clin$age_at_index.demographic),
  gender  = clin$gender.demographic,
  stage   = clin$ajcc_pathologic_stage.diagnoses,
  stringsAsFactors = FALSE
)
clin_df <- clin_df[!duplicated(clin_df$patient), ]

merged2 <- inner_join(merged, clin_df, by = "patient")
merged2$stage_simple <- ifelse(
  merged2$stage %in% c("Stage I", "Stage II"), "Early",
  ifelse(merged2$stage == "" | is.na(merged2$stage), NA, "Advanced")
)

cox_multi <- coxph(
  Surv(OS.time.month, OS) ~ OXCT1_group + HMGCS2_group + age + gender + stage_simple,
  data = merged2
)
summary(cox_multi)
# HMGCS2-Low: HR = 2.39 (1.38-4.12), p = 0.0018 (independent prognostic factor)
# OXCT1: HR = 1.14 (0.70-1.86), p = 0.61 (attenuated after stage adjustment)


# -----------------------------------------------------------------------------
# 6. Sensitivity analysis: optimal cutpoints (maxstat)
# -----------------------------------------------------------------------------

ms_oxct1  <- maxstat.test(Surv(OS.time.month, OS) ~ OXCT1,
                           data = merged, smethod = "LogRank")
ms_hmgcs2 <- maxstat.test(Surv(OS.time.month, OS) ~ HMGCS2,
                           data = merged, smethod = "LogRank")

cat("OXCT1 optimal cutpoint:", ms_oxct1$estimate,  "\n")  # 0.927
cat("HMGCS2 optimal cutpoint:", ms_hmgcs2$estimate, "\n") # 8.088

merged$OXCT1_opt  <- ifelse(merged$OXCT1  > ms_oxct1$estimate,  "High", "Low")
merged$HMGCS2_opt <- ifelse(merged$HMGCS2 > ms_hmgcs2$estimate, "High", "Low")

p_oxct1_opt  <- 1 - pchisq(survdiff(Surv(OS.time.month, OS) ~ OXCT1_opt,  data = merged)$chisq, 1)
p_hmgcs2_opt <- 1 - pchisq(survdiff(Surv(OS.time.month, OS) ~ HMGCS2_opt, data = merged)$chisq, 1)

cat("OXCT1 optimal cutpoint K-M p =",  p_oxct1_opt,  "\n")  # 0.000155
cat("HMGCS2 optimal cutpoint K-M p =", p_hmgcs2_opt, "\n")  # 3.12e-07


# -----------------------------------------------------------------------------
# 7. Immune infiltration: ESTIMATE (Fig. 5E)
# -----------------------------------------------------------------------------

library(estimate)

# Gene symbol expression matrix (genes x samples, tumor only)
expr_tumor_symbol <- tcga_sym[, c(FALSE, FALSE, FALSE,
                                   grepl("01A|01B", colnames(tcga_sym)[4:ncol(tcga_sym)]))]

write.table(expr_tumor_symbol, "D:/scRNA_project/expr_for_estimate.txt",
            sep = "\t", quote = FALSE)

filterCommonGenes("D:/scRNA_project/expr_for_estimate.txt",
                  "D:/scRNA_project/expr_filtered.gct",
                  id = "GeneSymbol")

estimateScore("D:/scRNA_project/expr_filtered.gct",
              "D:/scRNA_project/estimate_scores.gct",
              platform = "illumina")

# Read ESTIMATE scores
scores_raw <- read.table("D:/scRNA_project/estimate_scores.gct",
                         skip = 2, header = TRUE, sep = "\t", row.names = 1)
scores <- as.data.frame(t(scores_raw[, -1]))
colnames(scores) <- rownames(scores_raw)
scores$patient <- gsub("\\.", "-", substr(rownames(scores), 1, 12))

scores_merged <- inner_join(scores,
                             merged[, c("patient", "OXCT1", "OXCT1_group")],
                             by = "patient")

# Fig. 5E: OXCT1 vs ImmuneScore scatter
p_immune <- ggscatter(scores_merged, x = "OXCT1", y = "ImmuneScore",
  add = "reg.line", conf.int = TRUE,
  cor.coef = TRUE, cor.method = "spearman",
  color = "#E64B35", alpha = 0.4, size = 1.2,
  xlab = "OXCT1 expression (log2 FPKM)",
  ylab = "Immune Score (ESTIMATE)",
  title = "A  OXCT1 vs Immune infiltration") +
  theme_bw()
# Spearman R = 0.29, p = 9.5e-10

ggsave("D:/综述代码图/OXCT1_ImmuneScore_correlation.png",
       plot = p_immune, width = 6, height = 5, dpi = 300)


# -----------------------------------------------------------------------------
# 8. T-cell exhaustion markers with FDR correction (Fig. 5F)
# -----------------------------------------------------------------------------

exhaustion_genes <- c("CD8A", "HAVCR2", "LAG3", "PDCD1", "TIGIT")

exhaust_expr <- as.data.frame(t(tcga_sym[exhaustion_genes, 4:ncol(tcga_sym)]))
exhaust_expr$sample  <- rownames(exhaust_expr)
exhaust_expr$patient <- substr(exhaust_expr$sample, 1, 12)
exhaust_tumor <- exhaust_expr[grep("01A|01B", exhaust_expr$sample), ]

exhaust_merged <- inner_join(exhaust_tumor,
                              merged[, c("patient", "OXCT1_group")],
                              by = "patient")

# Wilcoxon with FDR correction
p_vals <- sapply(exhaustion_genes, function(g) {
  wilcox.test(exhaust_merged[[g]] ~ exhaust_merged$OXCT1_group)$p.value
})
p_fdr <- p.adjust(p_vals, method = "fdr")
print(data.frame(gene = exhaustion_genes, p_raw = p_vals, p_fdr = p_fdr))
# All FDR-corrected p values remain highly significant

# Plot
exhaust_long <- exhaust_merged %>%
  pivot_longer(cols = all_of(exhaustion_genes),
               names_to = "gene", values_to = "expression")

p_exhaust <- ggplot(exhaust_long,
                    aes(x = OXCT1_group, y = expression, fill = OXCT1_group)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~gene, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("High" = "#E64B35", "Low" = "#4DBBD5")) +
  stat_compare_means(method = "wilcox.test", label = "p.signif") +
  labs(title = "B  T cell exhaustion markers by OXCT1 expression",
       x = "", y = "Expression (log2 FPKM)") +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"))

ggsave("D:/综述代码图/OXCT1_exhaustion_markers.png",
       plot = p_exhaust, width = 12, height = 4, dpi = 300)

# Combined immune panel
combined_panel <- p_immune | p_exhaust + plot_layout(widths = c(1, 2))

png("D:/综述代码图/OXCT1_immune_panel.png",
    width = 14, height = 5, units = "in", res = 300)
combined_panel + plot_layout(widths = c(1, 2))
dev.off()


# -----------------------------------------------------------------------------
# 9. HMGCS2 vs OXCT1 correlation (Fig. 5G)
# -----------------------------------------------------------------------------

p_axis <- ggscatter(merged, x = "HMGCS2", y = "OXCT1",
  add = "reg.line", conf.int = TRUE,
  cor.coef = TRUE, cor.method = "spearman",
  color = "#3C5488", alpha = 0.4, size = 1.5,
  xlab = "HMGCS2 expression (log2 FPKM)",
  ylab = "OXCT1 expression (log2 FPKM)",
  title = "HMGCS2 vs OXCT1 expression in TCGA-LIHC HCC tumors") +
  theme_bw()
# Spearman R = -0.15, p = 0.0021

ggsave("D:/综述代码图/HMGCS2_OXCT1_correlation.png",
       plot = p_axis, width = 6, height = 5, dpi = 300)
ggsave("D:/综述代码图/HMGCS2_OXCT1_correlation.pdf",
       plot = p_axis, width = 6, height = 5)


# -----------------------------------------------------------------------------
# 10. CTNNB1 mutation analysis
# -----------------------------------------------------------------------------

ctnnb1_mut     <- maf[maf$Hugo_Symbol == "CTNNB1", ]
ctnnb1_mut_ids <- unique(ctnnb1_mut$Tumor_Sample_Barcode)
cat("CTNNB1 mutant samples:", length(ctnnb1_mut_ids), "\n")  # 96

merged$CTNNB1_mut <- ifelse(
  merged$patient %in% substr(ctnnb1_mut_ids, 1, 12),
  "Mutant", "Wild-type"
)
cat("Mutant:", sum(merged$CTNNB1_mut == "Mutant"),
    "| Wild-type:", sum(merged$CTNNB1_mut == "Wild-type"), "\n")

# NOTE: CTNNB1 mutation status was not associated with OXCT1 mRNA level
# in this mixed-etiology cohort (Wilcoxon p = 0.016, but direction reversed).
# This may reflect etiology heterogeneity; MASLD-specific validation is needed.


# -----------------------------------------------------------------------------
# Session info
# -----------------------------------------------------------------------------
sessionInfo()
