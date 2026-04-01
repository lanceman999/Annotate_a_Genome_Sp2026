library(dplyr)
library(data.table)

# Set working directory
setwd("/Users/lanceoconnor/Desktop/JohnsHopkins/classes/AnnotateAgenome/Annotate_a_Genome_Sp2026/lecture4_PRS")

# ============================================
# Load files
# ============================================

# Load cleaned discovery summary stats
sumstats <- fread("output/discovery_gwas_sim_sumstats_cleaned.tsv", header = TRUE)

cat("=== Discovery summary stats ===\n")
cat("Dimensions:", dim(sumstats), "\n")
head(sumstats)

# Load target BIM file
bim_path <- "/Users/lanceoconnor/Desktop/JohnsHopkins/classes/AnnotateAgenome/Annotate_a_Genome_Sp2026/exercise_4_GWAS_scripts/PopStrat_assignmentTwo/PLINK_files/QC/merged_1000G_HapMap.bim"

bim <- fread(bim_path, header = FALSE,
             col.names = c("CHR", "SNP", "CM", "BP", "A1_bim", "A2_bim"))

cat("\n=== Target BIM file ===\n")
cat("Dimensions:", dim(bim), "\n")
head(bim)

# ============================================
# 2a. SNP overlap with target
# ============================================

cat("\n=== 2a. SNP Overlap ===\n")

# Find overlapping SNPs
sumstats_snps <- sumstats$SNP
bim_snps <- bim$SNP

overlap_snps <- intersect(sumstats_snps, bim_snps)

cat("SNPs in discovery:", length(sumstats_snps), "\n")
cat("SNPs in target BIM:", length(bim_snps), "\n")
cat("Overlapping SNPs:", length(overlap_snps), "\n")
cat("% of discovery retained:", round(100 * length(overlap_snps) / length(sumstats_snps), 2), "%\n")
cat("% of target covered:", round(100 * length(overlap_snps) / length(bim_snps), 2), "%\n")

# Filter to overlapping SNPs
sumstats_overlap <- sumstats %>%
  filter(SNP %in% overlap_snps)

bim_overlap <- bim %>%
  filter(SNP %in% overlap_snps) %>%
  select(SNP, A1_bim, A2_bim)

# Merge sumstats with BIM alleles
harmonized <- sumstats_overlap %>%
  inner_join(bim_overlap, by = "SNP")

cat("\nRows after overlap filtering:", nrow(harmonized), "\n")

# ============================================
# 2b & 2c. Allele matching and flipping
# ============================================

cat("\n=== 2b & 2c. Allele Matching & Flipping ===\n")

# Define complement function for strand flipping
complement <- function(allele) {
  case_when(
    allele == "A" ~ "T",
    allele == "T" ~ "A",
    allele == "C" ~ "G",
    allele == "G" ~ "C",
    TRUE ~ NA_character_
  )
}

# Determine allele matching status
harmonized <- harmonized %>%
  mutate(
    # Original alleles from sumstats
    A1_ss = A1,
    A2_ss = A2,
    
    # Complement alleles
    A1_ss_comp = complement(A1),
    A2_ss_comp = complement(A2),
    
    # Check matching scenarios
    match_type = case_when(
      # Exact match (same order)
      (A1_ss == A1_bim & A2_ss == A2_bim) ~ "exact_match",
      
      # Swapped (A1/A2 reversed)
      (A1_ss == A2_bim & A2_ss == A1_bim) ~ "swapped",
      
      # Strand flip (complement, same order)
      (A1_ss_comp == A1_bim & A2_ss_comp == A2_bim) ~ "strand_flip",
      
      # Strand flip + swapped
      (A1_ss_comp == A2_bim & A2_ss_comp == A1_bim) ~ "strand_flip_swapped",
      
      # No match
      TRUE ~ "no_match"
    )
  )

# Report matching results
cat("\nAllele matching results:\n")
print(table(harmonized$match_type))

# ============================================
# 2d. Ambiguous strand SNPs (A/T and C/G)
# ============================================

cat("\n=== 2d. Ambiguous Strand SNPs ===\n")

harmonized <- harmonized %>%
  mutate(
    ambiguous = (A1_ss == "A" & A2_ss == "T") |
      (A1_ss == "T" & A2_ss == "A") |
      (A1_ss == "C" & A2_ss == "G") |
      (A1_ss == "G" & A2_ss == "C")
  )

n_ambiguous <- sum(harmonized$ambiguous)
cat("Ambiguous A/T and C/G SNPs:", n_ambiguous, "\n")

# Remove ambiguous SNPs
harmonized <- harmonized %>%
  filter(!ambiguous)

cat("Rows after removing ambiguous SNPs:", nrow(harmonized), "\n")

# ============================================
# Remove no-match SNPs
# ============================================

cat("\n=== Removing non-matching SNPs ===\n")

n_no_match <- sum(harmonized$match_type == "no_match")
cat("SNPs with no allele match:", n_no_match, "\n")

harmonized <- harmonized %>%
  filter(match_type != "no_match")

cat("Rows after removing no-match:", nrow(harmonized), "\n")

# ============================================
# Apply flipping logic
# ============================================

cat("\n=== Applying flip logic ===\n")

# Count SNPs that need flipping
n_swapped <- sum(harmonized$match_type %in% c("swapped", "strand_flip_swapped"))
cat("SNPs requiring BETA flip:", n_swapped, "\n")

# Apply harmonization
harmonized <- harmonized %>%
  mutate(
    # Flip BETA if alleles were swapped
    BETA_harmonized = case_when(
      match_type %in% c("swapped", "strand_flip_swapped") ~ -BETA,
      TRUE ~ BETA
    ),
    
    # Set A1 to match target BIM (effect allele aligned to target)
    A1_harmonized = A1_bim,
    A2_harmonized = A2_bim
  )

cat("\nBETA before vs after flipping (for swapped SNPs):\n")
harmonized %>%
  filter(match_type %in% c("swapped", "strand_flip_swapped")) %>%
  select(SNP, BETA, BETA_harmonized, match_type) %>%
  head(10) %>%
  print()

# ============================================
# Create final output file
# ============================================

cat("\n=== Creating output file ===\n")

prs_weights <- harmonized %>%
  select(
    SNP,
    A1 = A1_harmonized,    # Effect allele (aligned to target)
    BETA = BETA_harmonized, # Weight (sign flipped if needed)
    P
  )

cat("Final harmonized SNPs:", nrow(prs_weights), "\n")
head(prs_weights)

# Summary statistics
cat("\n=== Final summary ===\n")
cat("BETA range:", range(prs_weights$BETA), "\n")
cat("P range:", range(prs_weights$P), "\n")

# ============================================
# Save output
# ============================================

fwrite(prs_weights,
       file = "output/prs_weights_harmonized_qc.tsv",
       sep = "\t",
       quote = FALSE)

cat("\n✓ Saved: prs_weights_harmonized_qc.tsv\n")

# ============================================
# QC Summary Report
# ============================================

cat("\n")
cat("============================================\n")
cat("       HARMONIZATION SUMMARY REPORT        \n")
cat("============================================\n")
cat("Discovery SNPs (cleaned):", length(sumstats_snps), "\n")
cat("Target BIM SNPs:", length(bim_snps), "\n")
cat("Overlapping SNPs:", length(overlap_snps), "\n")
cat("Ambiguous SNPs removed:", n_ambiguous, "\n")
cat("Non-matching SNPs removed:", n_no_match, "\n")
cat("SNPs with flipped BETA:", n_swapped, "\n")
cat("Final harmonized SNPs:", nrow(prs_weights), "\n")
cat("% retained from discovery:", round(100 * nrow(prs_weights) / length(sumstats_snps), 2), "%\n")
cat("============================================\n")