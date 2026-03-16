library(ggplot2)
library(dplyr)

# ============================================
# Load MDS results
# ============================================
setwd("/Users/lanceoconnor/Desktop/JohnsHopkins/classes/AnnotateAgenome/Annotate_a_Genome_Sp2026/exercise_4_GWAS_scripts/PopStrat_assignmentTwo/PLINK_files/QC")

mds <- read.table("merged_mds.mds", header = TRUE)
head(mds)

# 1000 Genomes population info
panel_1000G <- read.table("/Users/lanceoconnor/Downloads/lecture_4/20100804.ALL.FIXED.panel", header = TRUE, stringsAsFactors = FALSE, fill = TRUE, col.names = c("sample", "pop", "platform1", "platform2")) %>%
  dplyr::select(sample,pop)
head(panel_1000G)

# Join by sample ID
mds_annotated <- mds %>%
  dplyr::left_join(panel_1000G, by = c("IID" = "sample"))  # adjust column name as needed

# Check for samples without population info (likely HapMap)
table(is.na(mds_annotated$pop))

# Basic plot: C1 vs C2
ggplot(mds_annotated, aes(x = C1, y = C2, color = pop)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(x = "MDS Dimension 1 (C1)",
       y = "MDS Dimension 2 (C2)",
       title = "MDS Analysis of 1000 Genomes and HapMap",
       color = "Population") +
  theme_bw() +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10))
