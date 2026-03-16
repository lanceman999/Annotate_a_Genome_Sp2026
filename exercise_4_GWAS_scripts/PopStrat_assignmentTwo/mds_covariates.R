
setwd("/Users/lanceoconnor/Desktop/JohnsHopkins/classes/AnnotateAgenome/Annotate_a_Genome_Sp2026/exercise_4_GWAS_scripts/PopStrat_assignmentTwo/PLINK_files/QC")

mds <- read.table("merged_mds.mds", header = TRUE)
head(mds)

# ============================================
# 7b) Create covariate file
# ============================================

# PLINK covariate format: FID, IID, C1, C2, ... C10
# The MDS file already has this format!

covariates <- mds[, c("FID", "IID", "C1", "C2", "C3", "C4", "C5", 
                      "C6", "C7", "C8", "C9", "C10")]

head(covariates)

write.table(covariates, file = "mds_covariates.txt",quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")