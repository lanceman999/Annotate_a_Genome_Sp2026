library(dplyr)
library(data.table)

setwd("/Users/lanceoconnor/Desktop/JohnsHopkins/classes/AnnotateAgenome/Annotate_a_Genome_Sp2026/lecture4_PRS")

# Load harmonized weights
prs_weights <- fread("output/prs_weights_harmonized_qc.tsv", header = TRUE)

# Load target BIM
bim_path <- "/Users/lanceoconnor/Desktop/JohnsHopkins/classes/AnnotateAgenome/Annotate_a_Genome_Sp2026/exercise_4_GWAS_scripts/PopStrat_assignmentTwo/PLINK_files/QC/merged_1000G_HapMap.bim"

bim <- fread(bim_path, header = FALSE,
             col.names = c("CHR", "SNP", "CM", "BP", "A1_bim", "A2_bim"))

# Load original cleaned sumstats (for comparison)
sumstats <- fread("output/discovery_gwas_sim_sumstats_cleaned.tsv", header = TRUE)

cat("============================================\n")
cat("    3) PRS-READY OUTPUT SANITY CHECKS      \n")
cat("============================================\n\n")

# ============================================
# 3a. Final weight file integrity
# ============================================

cat("=== 3a. Final Weight File Integrity ===\n\n")

# Check dimensions
cat("File dimensions:", dim(prs_weights), "\n")
cat("Columns:", paste(colnames(prs_weights), collapse = ", "), "\n\n")

# Check for unique SNP IDs
n_total <- nrow(prs_weights)
n_unique <- length(unique(prs_weights$SNP))
n_duplicates <- n_total - n_unique

cat("Total SNPs:", n_total, "\n")
cat("Unique SNPs:", n_unique, "\n")
cat("Duplicate SNPs:", n_duplicates, "\n")

if (n_duplicates == 0) {
  cat("✓ PASS: All SNP IDs are unique\n\n")
} else {
  cat("✗ FAIL: Duplicate SNP IDs found!\n")
  dup_snps <- prs_weights$SNP[duplicated(prs_weights$SNP)]
  cat("Duplicates:", head(dup_snps, 10), "\n\n")
}

# Check for missing BETA
n_missing_beta <- sum(is.na(prs_weights$BETA))
cat("Missing BETA values:", n_missing_beta, "\n")

if (n_missing_beta == 0) {
  cat("✓ PASS: No missing BETA values\n\n")
} else {
  cat("✗ FAIL: Missing BETA values found!\n\n")
}

# Check for missing A1
n_missing_a1 <- sum(is.na(prs_weights$A1))
cat("Missing A1 values:", n_missing_a1, "\n")

if (n_missing_a1 == 0) {
  cat("✓ PASS: No missing A1 values\n\n")
} else {
  cat("✗ FAIL: Missing A1 values found!\n\n")
}

# Check for missing P
n_missing_p <- sum(is.na(prs_weights$P))
cat("Missing P values:", n_missing_p, "\n")

if (n_missing_p == 0) {
  cat("✓ PASS: No missing P values\n\n")
} else {
  cat("✗ FAIL: Missing P values found!\n\n")
}

# BETA summary statistics
cat("BETA summary:\n")
print(summary(prs_weights$BETA))
cat("\n")

# P-value summary
cat("P-value summary:\n")
print(summary(prs_weights$P))
cat("\n")

# ============================================
# 3b. Spot-check allele alignment
# ============================================

cat("=== 3b. Spot-Check Allele Alignment ===\n\n")

# Randomly sample 10 SNPs
set.seed(42)  # For reproducibility
sample_snps <- sample(prs_weights$SNP, 10)

# Get info from all sources
spot_check <- prs_weights %>%
  filter(SNP %in% sample_snps) %>%
  rename(A1_harmonized = A1, BETA_harmonized = BETA) %>%
  left_join(
    bim %>% select(SNP, A1_bim, A2_bim),
    by = "SNP"
  ) %>%
  left_join(
    sumstats %>% select(SNP, A1, A2, BETA) %>% 
      rename(A1_original = A1, A2_original = A2, BETA_original = BETA),
    by = "SNP"
  ) %>%
  mutate(
    # Check if BETA was flipped
    beta_flipped = (BETA_harmonized != BETA_original),
    
    # Check if A1_harmonized matches BIM
    a1_matches_bim = (A1_harmonized == A1_bim)
  ) %>%
  select(SNP, 
         A1_original, A2_original, BETA_original,
         A1_bim, A2_bim,
         A1_harmonized, BETA_harmonized,
         beta_flipped, a1_matches_bim)

cat("Randomly sampled 10 SNPs for manual verification:\n\n")
print(as.data.frame(spot_check))

cat("\n")
cat("Column descriptions:\n")
cat("  A1_original, A2_original: Alleles from discovery sumstats\n")
cat("  BETA_original: Original effect size from discovery\n")
cat("  A1_bim, A2_bim: Alleles from target BIM file\n")
cat("  A1_harmonized: Effect allele in output (aligned to target)\n")
cat("  BETA_harmonized: Effect size in output (flipped if needed)\n")
cat("  beta_flipped: TRUE if sign was flipped\n")
cat("  a1_matches_bim: TRUE if output A1 matches BIM A1\n")
cat("\n")

# Verify alignment logic
cat("=== Verification Logic ===\n\n")

for (i in 1:nrow(spot_check)) {
  snp <- spot_check$SNP[i]
  cat("SNP:", snp, "\n")
  cat("  Original: A1=", spot_check$A1_original[i], 
      ", A2=", spot_check$A2_original[i],
      ", BETA=", round(spot_check$BETA_original[i], 4), "\n", sep = "")
  cat("  BIM:      A1=", spot_check$A1_bim[i],
      ", A2=", spot_check$A2_bim[i], "\n", sep = "")
  cat("  Output:   A1=", spot_check$A1_harmonized[i],
      ", BETA=", round(spot_check$BETA_harmonized[i], 4), "\n", sep = "")
  
  # Check logic
  if (spot_check$A1_original[i] == spot_check$A1_bim[i]) {
    expected_flip <- FALSE
    cat("  Expected: No flip (alleles match)\n")
  } else if (spot_check$A1_original[i] == spot_check$A2_bim[i]) {
    expected_flip <- TRUE
    cat("  Expected: Flip (alleles swapped)\n")
  } else {
    cat("  Expected: Strand flip scenario\n")
    expected_flip <- NA
  }
  
  if (!is.na(expected_flip)) {
    if (spot_check$beta_flipped[i] == expected_flip) {
      cat("  ✓ CORRECT\n")
    } else {
      cat("  ✗ CHECK NEEDED\n")
    }
  }
  cat("\n")
}

# ============================================
# Summary
# ============================================

cat("============================================\n")
cat("           SANITY CHECK SUMMARY            \n")
cat("============================================\n")
cat("Total SNPs in output:", nrow(prs_weights), "\n")
cat("Unique SNPs: ✓" , ifelse(n_duplicates == 0, "PASS", "FAIL"), "\n")
cat("No missing BETA:", ifelse(n_missing_beta == 0, "✓ PASS", "✗ FAIL"), "\n")
cat("No missing A1:", ifelse(n_missing_a1 == 0, "✓ PASS", "✗ FAIL"), "\n")
cat("No missing P:", ifelse(n_missing_p == 0, "✓ PASS", "✗ FAIL"), "\n")
cat("Spot-check complete: 10 SNPs manually verified\n")
cat("============================================\n")

# ============================================
# Save QC'd file (same as before, confirms it's ready)
# ============================================

# The file is already QC'd, but we can rename/confirm
file.copy("output/prs_weights_harmonized_qc.tsv", "output/prs_weights_harmonized_qcd.tsv", overwrite = TRUE)
cat("\n✓ Confirmed: prs_weights_harmonized_qcd.tsv\n")