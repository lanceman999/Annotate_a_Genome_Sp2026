library(dplyr)
library(data.table)

setwd("/Users/lanceoconnor/Desktop/JohnsHopkins/classes/AnnotateAgenome/Annotate_a_Genome_Sp2026/lecture4_PRS")

# ============================================
# 1. Load discovery GWAS summary stats
# ============================================

sumstats <- fread("discovery_gwas_sim_sumstats.tsv.gz", header = TRUE)

cat("=== Original file ===\n")
cat("Dimensions:", dim(sumstats), "\n")
head(sumstats)
colnames(sumstats)

# ============================================
# 1a. Column/format validation
# ============================================

required_cols <- c("SNP", "A1", "A2", "BETA", "SE", "P", "N", "CHR", "BP")

cat("\n=== Column validation ===\n")
cat("Required columns:", paste(required_cols, collapse = ", "), "\n")
cat("Present columns:", paste(colnames(sumstats), collapse = ", "), "\n")

missing_cols <- setdiff(required_cols, colnames(sumstats))
if (length(missing_cols) == 0) {
  cat("✓ All required columns present!\n")
} else {
  cat("✗ Missing columns:", paste(missing_cols, collapse = ", "), "\n")
}

# ============================================
# 1b. Missingness check
# ============================================

cat("\n=== Missingness check ===\n")

# Check missingness for each required column
for (col in required_cols) {
  if (col %in% colnames(sumstats)) {
    n_missing <- sum(is.na(sumstats[[col]]))
    pct_missing <- round(100 * n_missing / nrow(sumstats), 3)
    cat(col, ": ", n_missing, " missing (", pct_missing, "%)\n", sep = "")
  }
}

# Total rows with any missing values in key columns
rows_before <- nrow(sumstats)

sumstats_clean <- sumstats %>%
  filter(!is.na(SNP),
         !is.na(A1),
         !is.na(A2),
         !is.na(BETA),
         !is.na(SE),
         !is.na(P),
         !is.na(N))

rows_after_missing <- nrow(sumstats_clean)
cat("\nRows removed due to missingness:", rows_before - rows_after_missing, "\n")

# ============================================
# 1c. Range checks
# ============================================

cat("\n=== Range checks ===\n")

# Check P-value range (0 ≤ P ≤ 1)
p_out_of_range <- sum(sumstats_clean$P < 0 | sumstats_clean$P > 1)
cat("P-values out of range [0,1]:", p_out_of_range, "\n")

# Remove invalid P-values
sumstats_clean <- sumstats_clean %>%
  filter(P >= 0 & P <= 1)

# Check SE > 0
se_invalid <- sum(sumstats_clean$SE <= 0)
cat("SE values ≤ 0:", se_invalid, "\n")

# Remove invalid SE
sumstats_clean <- sumstats_clean %>%
  filter(SE > 0)

# Flag extreme BETA outliers
beta_mean <- mean(sumstats_clean$BETA)
beta_sd <- sd(sumstats_clean$BETA)
beta_threshold <- 5 * beta_sd  # Flag if |BETA| > 5 SD from mean

extreme_beta <- sum(abs(sumstats_clean$BETA - beta_mean) > beta_threshold)
cat("Extreme |BETA| outliers (>5 SD):", extreme_beta, "\n")

# View extreme outliers
if (extreme_beta > 0) {
  cat("\nExtreme BETA values:\n")
  sumstats_clean %>%
    filter(abs(BETA - beta_mean) > beta_threshold) %>%
    select(SNP, CHR, BP, BETA, SE, P) %>%
    print()
}


rows_after_range <- nrow(sumstats_clean)
cat("\nRows after range checks:", rows_after_range, "\n")

# ============================================
# 1d. Duplicate SNP IDs
# ============================================

cat("\n=== Duplicate SNP check ===\n")

# Count duplicates
n_duplicates <- sum(duplicated(sumstats_clean$SNP))
cat("Duplicate SNP IDs:", n_duplicates, "\n")

# Remove duplicates, keeping the one with smallest P-value
sumstats_clean <- sumstats_clean %>%
  arrange(P) %>%
  distinct(SNP, .keep_all = TRUE)

rows_after_dedup <- nrow(sumstats_clean)
cat("Rows after removing duplicates:", rows_after_dedup, "\n")

# ============================================
# Summary
# ============================================

cat("\n=== QC Summary ===\n")
cat("Original rows:", rows_before, "\n")
cat("Removed - missingness:", rows_before - rows_after_missing, "\n")
cat("Removed - range checks:", rows_after_missing - rows_after_range, "\n")
cat("Removed - duplicates:", rows_after_range - rows_after_dedup, "\n")
cat("Final rows:", rows_after_dedup, "\n")

# ============================================
# Save cleaned file
# ============================================

fwrite(sumstats_clean, 
       file = "discovery_gwas_sim_sumstats_cleaned.tsv",
       sep = "\t",
       quote = FALSE)

cat("\n✓ Saved: discovery_gwas_sim_sumstats_cleaned.tsv\n")

# ============================================
# Final validation
# ============================================

cat("\n=== Final file check ===\n")
cat("Dimensions:", dim(sumstats_clean), "\n")
head(sumstats_clean)
summary(sumstats_clean$P)
summary(sumstats_clean$BETA)
summary(sumstats_clean$SE)