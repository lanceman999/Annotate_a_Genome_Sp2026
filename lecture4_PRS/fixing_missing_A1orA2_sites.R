library(data.table)
library(dplyr)

setwd("/Users/lanceoconnor/Desktop/JohnsHopkins/classes/AnnotateAgenome/Annotate_a_Genome_Sp2026/lecture4_PRS/output")

# Load file
sumstats <- fread("discovery_gwas_sim_sumstats_cleaned.tsv", header = TRUE)

cat("=== Checking for missing/empty alleles ===\n")
cat("Original rows:", nrow(sumstats), "\n\n")

# Check for empty or NA alleles
missing_a1 <- sum(is.na(sumstats$A1) | sumstats$A1 == "")
missing_a2 <- sum(is.na(sumstats$A2) | sumstats$A2 == "")

cat("Missing/empty A1:", missing_a1, "\n")
cat("Missing/empty A2:", missing_a2, "\n")

# Show problematic rows
cat("\nProblematic rows:\n")
problematic <- sumstats %>%
  filter(is.na(A1) | A1 == "" | is.na(A2) | A2 == "")
print(problematic)

# ============================================
# Remove rows with missing alleles
# ============================================

sumstats_fixed <- sumstats %>%
  filter(!is.na(A1) & A1 != "",
         !is.na(A2) & A2 != "",
         !is.na(P) & P >= 0 & P <= 1,
         !is.na(BETA),
         !is.na(SE) & SE > 0)

cat("\n=== After fixing ===\n")
cat("Rows remaining:", nrow(sumstats_fixed), "\n")
cat("Rows removed:", nrow(sumstats) - nrow(sumstats_fixed), "\n")

# Verify P-value range
cat("\nP-value range:", range(sumstats_fixed$P), "\n")

# ============================================
# Save fixed file
# ============================================

fwrite(sumstats_fixed, 
       file = "discovery_gwas_sim_sumstats_cleaned.tsv",
       sep = "\t",
       quote = FALSE)

cat("\n✓ Saved: discovery_gwas_sim_sumstats_cleaned.tsv\n")