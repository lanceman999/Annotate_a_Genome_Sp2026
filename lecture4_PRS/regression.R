# ============================================
# PRS Assignment - Question 4.4
# Logistic Regression: phenotype ~ covariates + PRS
# ============================================

library(data.table)
library(dplyr)

# Set working directory
setwd("/Users/lanceoconnor/Desktop/JohnsHopkins/classes/AnnotateAgenome/Annotate_a_Genome_Sp2026/lecture4_PRS/output")

# ============================================
# Load phenotype and covariate files
# ============================================

cat("=== Loading files ===\n\n")

# Load phenotype
pheno <- fread("pheno.tsv", header = TRUE) %>%
  mutate(FID = as.character(FID),
         IID = as.character(IID))

cat("Phenotype file:", nrow(pheno), "individuals\n")
cat("Phenotype distribution:\n")
print(table(pheno$PHENO))

# Load covariates (MDS components)
covar <- fread("covar.tsv", header = TRUE) %>%
  mutate(FID = as.character(FID),
         IID = as.character(IID))

cat("\nCovariate file:", nrow(covar), "individuals\n")
cat("Covariates:", paste(colnames(covar)[3:ncol(covar)], collapse = ", "), "\n")

# ============================================
# Define thresholds
# ============================================

thresholds <- c("5e-8", "1e-6", "1e-4", "1e-3", "1e-2", "0.05", "0.1", "0.5", "1")

# ============================================
# For each threshold, run logistic regression
# ============================================

cat("\n=== Running Logistic Regression for Each Threshold ===\n\n")

# Store results
results_summary <- data.frame(
  threshold = character(),
  n_samples = integer(),
  n_snps = integer(),
  prs_beta = numeric(),
  prs_se = numeric(),
  prs_z = numeric(),
  prs_pvalue = numeric(),
  stringsAsFactors = FALSE
)

for (thresh in thresholds) {
  
  cat("----------------------------------------\n")
  cat("Threshold: P ≤", thresh, "\n")
  cat("----------------------------------------\n")
  
  # ============================================
  # Step 1: Load PRS profile
  # ============================================
  
  profile_file <- paste0("prs_", thresh, ".profile")
  
  if (!file.exists(profile_file)) {
    cat("Warning: ", profile_file, " not found, skipping...\n\n")
    next
  }
  
  prs <- fread(profile_file, header = TRUE)
  cat("PRS file:", nrow(prs), "individuals\n")
  
  # Get number of SNPs used
  snp_file <- paste0("prs_", thresh, ".snps")
  n_snps <- length(readLines(snp_file))
  cat("SNPs used:", n_snps, "\n")
  
  # ============================================
  # Step 2: Merge phenotype + covariates + PRS
  # ============================================
  
  # Select relevant columns from PRS and convert types
  prs_slim <- prs %>%
    mutate(FID = as.character(FID),
           IID = as.character(IID)) %>%
    select(FID, IID, SCORE) %>%
    rename(PRS = SCORE)
  
  # Merge all data
  merged <- pheno %>%
    inner_join(covar, by = c("FID", "IID")) %>%
    inner_join(prs_slim, by = c("FID", "IID"))
  
  # Remove missing phenotypes (-9)
  merged <- merged %>%
    filter(PHENO != -9)
  
  # Convert phenotype to 0/1 for logistic regression
  # PLINK coding: 1 = control, 2 = case
  # R glm coding: 0 = control, 1 = case
  merged <- merged %>%
    mutate(PHENO_binary = PHENO - 1)
  
  cat("Merged samples (with phenotype):", nrow(merged), "\n")
  cat("Cases:", sum(merged$PHENO_binary == 1), "\n")
  cat("Controls:", sum(merged$PHENO_binary == 0), "\n")
  
  # ============================================
  # Step 3: Standardize PRS
  # ============================================
  
  merged <- merged %>%
    mutate(PRS_scaled = scale(PRS)[,1])
  
  # ============================================
  # Step 4: Run logistic regression
  # ============================================
  
  # Get covariate column names
  covar_cols <- colnames(covar)[3:ncol(covar)]
  
  formula_str <- paste("PHENO_binary ~ PRS_scaled +", paste(covar_cols, collapse = " + "))
  cat("Model:", formula_str, "\n")
  
  model <- glm(as.formula(formula_str), 
               data = merged, 
               family = binomial(link = "logit"))
  
  # ============================================
  # Step 5: Extract results
  # ============================================
  
  model_summary <- summary(model)
  coef_table <- as.data.frame(model_summary$coefficients)
  
  # Get PRS coefficient
  prs_coef <- coef_table["PRS_scaled", ]
  
  cat("\nPRS Association Results:\n")
  cat("  Beta (log-OR):", round(prs_coef$Estimate, 4), "\n")
  cat("  SE:", round(prs_coef$`Std. Error`, 4), "\n")
  cat("  Z:", round(prs_coef$`z value`, 4), "\n")
  cat("  P-value:", format(prs_coef$`Pr(>|z|)`, digits = 4), "\n")
  
  # Odds ratio
  or <- exp(prs_coef$Estimate)
  or_ci_low <- exp(prs_coef$Estimate - 1.96 * prs_coef$`Std. Error`)
  or_ci_high <- exp(prs_coef$Estimate + 1.96 * prs_coef$`Std. Error`)
  cat("  OR (95% CI):", round(or, 3), "(", round(or_ci_low, 3), "-", round(or_ci_high, 3), ")\n")
  
  # ============================================
  # Step 6: Save full results to file
  # ============================================
  
  output_file <- paste0("assoc_$", thresh, ".assoc.logistic")
  
  # Create output table
  output_table <- data.frame(
    VARIABLE = rownames(coef_table),
    BETA = coef_table$Estimate,
    SE = coef_table$`Std. Error`,
    Z = coef_table$`z value`,
    P = coef_table$`Pr(>|z|)`
  )
  
  fwrite(output_table, file = output_file, sep = "\t", quote = FALSE)
  cat("\nSaved:", output_file, "\n\n")
  
  # Store summary results
  results_summary <- rbind(results_summary, data.frame(
    threshold = thresh,
    n_samples = nrow(merged),
    n_snps = n_snps,
    prs_beta = prs_coef$Estimate,
    prs_se = prs_coef$`Std. Error`,
    prs_z = prs_coef$`z value`,
    prs_pvalue = prs_coef$`Pr(>|z|)`,
    stringsAsFactors = FALSE
  ))
}

# ============================================
# Summary across all thresholds
# ============================================

cat("\n")
cat("============================================\n")
cat("     SUMMARY ACROSS ALL THRESHOLDS         \n")
cat("============================================\n\n")

results_summary <- results_summary %>%
  mutate(
    OR = exp(prs_beta),
    significant = ifelse(prs_pvalue < 0.05, "*", "")
  )

print(results_summary %>% 
        select(threshold, n_snps, prs_beta, prs_se, prs_pvalue, OR, significant) %>%
        mutate(across(c(prs_beta, prs_se, OR), ~round(.x, 4)),
               prs_pvalue = format(prs_pvalue, digits = 3)))

# ============================================
# Find best threshold
# ============================================

best_thresh <- results_summary %>%
  filter(prs_pvalue == min(prs_pvalue))

cat("\n")
cat("Best threshold (lowest P-value):\n")
cat("  Threshold:", best_thresh$threshold, "\n")
cat("  N SNPs:", best_thresh$n_snps, "\n")
cat("  P-value:", format(best_thresh$prs_pvalue, digits = 4), "\n")
cat("  OR:", round(best_thresh$OR, 3), "\n")

# ============================================
# Save summary table
# ============================================

fwrite(results_summary, file = "prs_association_summary.tsv", sep = "\t", quote = FALSE)
cat("\n✓ Saved: prs_association_summary.tsv\n")

# ============================================
# Verify output files
# ============================================

cat("\n=== Output files created ===\n")
assoc_files <- list.files(pattern = "^assoc_.*\\.assoc\\.logistic$")
for (f in assoc_files) {
  cat(" ", f, "\n")
}

cat("\n=== Done! ===\n")