library(dplyr)

setwd("/Users/lanceoconnor/Desktop/JohnsHopkins/classes/AnnotateAgenome/Annotate_a_Genome_Sp2026/exercise_4_GWAS_scripts/GWAS_assignmentThree/output")

logistic <- read.table("logistic_assoc.assoc.logistic", header = TRUE, stringsAsFactors = FALSE)

head(logistic)

# ============================================
# Filter to keep only the SNP results (TEST == "ADD")
# This removes the covariate rows (C1, C2, etc.)
# ============================================

logistic_snps <- logistic %>%
  dplyr::filter(TEST == "ADD")

head(logistic_snps)
nrow(logistic_snps) # 1102192

# ============================================
# Remove rows with NA p-values
# ============================================

logistic_clean <- logistic_snps %>%
  dplyr::filter(!is.na(P))

nrow(logistic_clean)

write.table(logistic_clean, file = "logistic_clean.txt", quote = FALSE, row.names = FALSE, sep = "\t")