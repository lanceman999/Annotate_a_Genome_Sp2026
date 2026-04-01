library(data.table)
library(dplyr)

setwd("/Users/lanceoconnor/Desktop/JohnsHopkins/classes/AnnotateAgenome/Annotate_a_Genome_Sp2026/lecture4_PRS/output")

# Load cleaned summary stats
sumstats <- fread("discovery_gwas_sim_sumstats_cleaned.tsv", header = TRUE)
cat("Loaded summary stats:", nrow(sumstats), "SNPs\n")

# Load clumped SNPs
clumped_snps <- fread("clumped.snps", header = FALSE)$V1
cat("Loaded clumped SNPs:", length(clumped_snps), "SNPs\n\n")

# ============================================
# Define P-value thresholds
# ============================================

thresholds <- c(5e-8, 1e-6, 1e-4, 1e-3, 1e-2, 0.05, 0.1, 0.5, 1.0)

cat("P-value thresholds:", paste(thresholds, collapse = ", "), "\n\n")

# ============================================
# For each threshold, create SNP set
# ============================================

cat("=== Creating thresholded SNP sets ===\n\n")

# Store results for summary
results <- data.frame(
  threshold = thresholds,
  n_pass_pvalue = NA,
  n_clumped = NA,
  n_final = NA
)

for (i in seq_along(thresholds)) {
  
  thresh <- thresholds[i]
  
  # Step 1: SNPs passing P-value threshold
  snps_pass_p <- sumstats %>%
    filter(P <= thresh) %>%
    pull(SNP)
  
  # Step 2: Intersect with clumped SNPs
  snps_final <- intersect(snps_pass_p, clumped_snps)
  
  # Step 3: Save to file
  # Format threshold for filename
  if (thresh < 0.01) {
    thresh_name <- format(thresh, scientific = TRUE)
    thresh_name <- gsub("e-0", "e-", thresh_name)  # 5e-08 -> 5e-8
    thresh_name <- gsub("e\\+", "e", thresh_name)
  } else {
    thresh_name <- as.character(thresh)
  }
  
  filename <- paste0("prs_", thresh_name, ".snps")
  
  # Write SNP list
  writeLines(snps_final, filename)
  
  # Store results
  results$n_pass_pvalue[i] <- length(snps_pass_p)
  results$n_clumped[i] <- length(clumped_snps)
  results$n_final[i] <- length(snps_final)
  
  cat(sprintf("Threshold P ≤ %s:\n", thresh_name))
  cat(sprintf("  SNPs passing P-value: %d\n", length(snps_pass_p)))
  cat(sprintf("  After clumping intersection: %d\n", length(snps_final)))
  cat(sprintf("  Saved to: %s\n\n", filename))
}

# ============================================
# Summary table
# ============================================

cat("=== Summary Table ===\n\n")
print(results)

# ============================================
# Verify output files
# ============================================

cat("\n=== Output files created ===\n")
snp_files <- list.files(pattern = "^prs_.*\\.snps$")
for (f in snp_files) {
  n_snps <- length(readLines(f))
  cat(sprintf("  %s: %d SNPs\n", f, n_snps))
}

cat("\n✓ Done! All prs_{T}.snps files created.\n")