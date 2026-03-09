library(ggplot2)
library(readr)
library(dplyr)
library(cowplot)

###################################################################################################
###################################################################################################

# SNP MISSINGNESS #

###################################################################################################
###################################################################################################
# Individual SNP missingness
imiss <- read_table("/Users/lanceoconnor/Desktop/JohnsHopkins/classes/AnnotateAgenome/Annotate_a_Genome_Sp2026/exercise_4_GWAS_scripts/missing_snps/missing_SNPs_out.imiss")

ggplot(imiss, aes(x = F_MISS)) +
  geom_histogram(binwidth = 0.005, fill = "blue", color = "black") +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text = element_text(size = 12, color ='black'),
    axis.title = element_text(size = 14, color = 'black'),
    plot.title = element_text(size = 18, color = 'black', hjust = 0.5)
  ) +
  scale_y_log10() +
  labs(x = "Individual Missingness (F_MISS)", y = "Number of Individuals", 
       title = "Histogram of Individual Missingness") +
  scale_y_continuous(expand = c(0,2)) +
  scale_x_continuous(expan = c(0,0.0005))


# SNP missingness
lmiss <- read_table("/Users/lanceoconnor/Desktop/JohnsHopkins/classes/AnnotateAgenome/Annotate_a_Genome_Sp2026/exercise_4_GWAS_scripts/missing_snps/missing_SNPs_out.lmiss")

ggplot(lmiss, aes(x = F_MISS)) +
  geom_histogram(binwidth = 0.005, fill = "firebrick", color = "black") +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text = element_text(size = 12, color ='black'),
    axis.title = element_text(size = 14, color = 'black'),
    plot.title = element_text(size = 18, color = 'black', hjust = 0.5)
  ) +
  labs(x = "SNP Missingness (F_MISS)", y = "Number of SNPs", title = "Histogram of SNP Missingness") +
  scale_y_continuous(expand = c(0,10000)) +
  scale_x_continuous(expan = c(0,0.0005))





###################################################################################################
###################################################################################################

# SEX DISCREPANCIES #

###################################################################################################
###################################################################################################
sex_check <- readr::read_tsv("/Users/lanceoconnor/Desktop/JohnsHopkins/classes/AnnotateAgenome/Annotate_a_Genome_Sp2026/exercise_4_GWAS_scripts/sex_discrep/sex_check.tabbed.tsv")

# Read FAM file
fam <- read.table('/Users/lanceoconnor/Downloads/lecture_4/HapMap_3_r3_1.fam', header = FALSE, 
                  col.names = c('FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENO'))

# Merge
final_plt <- sex_check %>%
  dplyr::left_join(fam[, c('FID', 'IID', 'SEX')], by = c('FID', 'IID')) %>%
  dplyr::mutate(sex_label = ifelse(SEX == 1, 'Male', 'Female'))

# Create figure with 3 plots
p1 <- ggplot(data = final_plt, aes(x = F)) +
  geom_histogram(bins = 50, fill = 'blue', alpha = 0.7, color = 'black') +
  labs(title = paste('All Samples (N =', nrow(final_plt), ')'),
       y = 'Frequency') +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'),
        axis.title.x = element_blank())

p2 <- ggplot(dplyr::filter(final_plt, SEX == 1), aes(x = F)) +
  geom_histogram(bins = 50, fill = 'green', alpha = 0.7, color = 'black') +
  labs(title = paste('Reported Males (N =', nrow(dplyr::filter(final_plt, SEX == 1)), ')'),
       y = 'Frequency') +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'),
        axis.title.x = element_blank()) 

p3 <- ggplot(dplyr::filter(final_plt, SEX == 2), aes(x = F)) +
  geom_histogram(bins = 50, fill = 'red', alpha = 0.7, color = 'black') +
  labs(title = paste('Reported Females (N =', nrow(dplyr::filter(final_plt, SEX == 2)), ')'),
       x = 'Inbreeding Coefficient (F)',
       y = 'Frequency') +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))

concatenated_hist <- cowplot::plot_grid(
  p1,p2,p3,
  nrow = 3
)
concatenated_hist






###################################################################################################
###################################################################################################

# ALLELE FREQUENCY #

###################################################################################################
###################################################################################################
af <- read.table("/Users/lanceoconnor/Desktop/JohnsHopkins/classes/AnnotateAgenome/Annotate_a_Genome_Sp2026/exercise_4_GWAS_scripts/MAF/allele_freq.frq", header = TRUE)

ggplot(af, aes(x = MAF)) +
  geom_histogram(bins = 100, fill = "steelblue", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 0.01, color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "Minor Allele Frequency Distribution",
       x = "MAF",
       y = "Count",
       caption = "Red line = MAF threshold (0.01)") +
  theme_bw() + 
  theme( panel.grid.major.x = element_blank(),
         panel.grid.minor.x = element_blank(),
         axis.text = element_text(size = 12, color ='black'),
         axis.title = element_text(size = 14, color = 'black'),
         plot.title = element_text(size = 18, color = 'black', hjust = 0.5)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))


ggplot(af %>% dplyr::filter(MAF > 0.01), aes(x = MAF)) +
  geom_histogram(bins = 100, fill = "steelblue", color = "black", alpha = 0.7) +
  # geom_vline(xintercept = 0.01, color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "Minor Allele Frequency Distribution",
       x = "MAF",
       y = "Count") +
  theme_bw() + 
  theme( panel.grid.major.x = element_blank(),
         panel.grid.minor.x = element_blank(),
         axis.text = element_text(size = 12, color ='black'),
         axis.title = element_text(size = 14, color = 'black'),
         plot.title = element_text(size = 18, color = 'black', hjust = 0.5)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0))







###################################################################################################
###################################################################################################

# HARDY-WEINBERG #

###################################################################################################
###################################################################################################
hwe <- read.table("/Users/lanceoconnor/Desktop/JohnsHopkins/classes/AnnotateAgenome/Annotate_a_Genome_Sp2026/exercise_4_GWAS_scripts/HWE/HWE_filtered.hwe", head = TRUE) %>%
  dplyr::filter(TEST == "ALL")


ggplot(hwe, aes(x = -log10(P))) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black", alpha = 0.7) +
  geom_vline(xintercept = -log10(1e-6), color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "HWE -log10(P-value) Distribution (After Filtering)",
       x = "-log10(HWE P-value)",
       y = "Count",
       caption = "Red line = HWE threshold (1e-6)") +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text = element_text(size = 12, color ='black'),
    axis.title = element_text(size = 14, color = 'black'),
    plot.title = element_text(size = 18, color = 'black', hjust = 0.5))






###################################################################################################
###################################################################################################

# HETEROZYGOSITY #

###################################################################################################
###################################################################################################
het <- read.table('/Users/lanceoconnor/Desktop/JohnsHopkins/classes/AnnotateAgenome/Annotate_a_Genome_Sp2026/exercise_4_GWAS_scripts/heterozygosity/het_stats.het', header = TRUE)

het <- het %>%
  dplyr::mutate(het_rate = (N.NM. - O.HOM.) / N.NM.)

het_mean <- mean(het$het_rate)
het_sd <- sd(het$het_rate)

het <- het %>%
  dplyr::mutate(
    lower_bound = het_mean - 3 * het_sd,
    upper_bound = het_mean + 3 * het_sd,
    outlier = het_rate < lower_bound | het_rate > upper_bound)

outliers <- het %>% dplyr::filter(outlier == TRUE)

cat("=== HETEROZYGOSITY OUTLIERS ===\n")
cat("Mean het rate:", het_mean, "\n")
cat("SD het rate:", het_sd, "\n")
cat("Lower bound (mean - 3SD):", het_mean - 3 * het_sd, "\n")
cat("Upper bound (mean + 3SD):", het_mean + 3 * het_sd, "\n")
cat("Number of outliers:", nrow(outliers), "\n\n")

if (nrow(outliers) > 0) {
  cat("Outlier individuals:\n")
  print(outliers %>% select(FID, IID, het_rate, F))
}

write.table(outliers %>% dplyr::select(FID, IID), "/Users/lanceoconnor/Desktop/JohnsHopkins/classes/AnnotateAgenome/Annotate_a_Genome_Sp2026/exercise_4_GWAS_scripts/heterozygosity/het_outliers.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Plot histogram
ggplot(het, aes(x = het_rate)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black", alpha = 0.7) +
  geom_vline(xintercept = het_mean, color = "black", linetype = "solid", linewidth = 1) +
  geom_vline(xintercept = het_mean - 3 * het_sd, color = "red", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = het_mean + 3 * het_sd, color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "Heterozygosity Rate Distribution",
       x = "Heterozygosity Rate",
       y = "Count",
       caption = "Black line = mean; Red lines = ±3 SD threshold") +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text = element_text(size = 12, color ='black'),
    axis.title = element_text(size = 14, color = 'black'),
    plot.title = element_text(size = 18, color = 'black', hjust = 0.5))

######### After filtering
het_filtered <- read.table("/Users/lanceoconnor/Desktop/JohnsHopkins/classes/AnnotateAgenome/Annotate_a_Genome_Sp2026/exercise_4_GWAS_scripts/heterozygosity/het_stats_filtered.het", header = TRUE)

het_filtered <- het_filtered %>%
  mutate(het_rate = (N.NM. - O.HOM.) / N.NM.)

het_mean <- mean(het_filtered$het_rate)
het_sd <- sd(het_filtered$het_rate)

ggplot(het_filtered, aes(x = het_rate)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black", alpha = 0.7) +
  geom_vline(xintercept = het_mean, color = "black", linetype = "solid", linewidth = 1) +
  geom_vline(xintercept = het_mean - 3 * het_sd, color = "red", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = het_mean + 3 * het_sd, color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "Heterozygosity Rate Distribution (After Filtering)",
       x = "Heterozygosity Rate",
       y = "Count",
       caption = "Black line = mean; Red lines = ±3 SD threshold") +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text = element_text(size = 12, color ='black'),
    axis.title = element_text(size = 14, color = 'black'),
    plot.title = element_text(size = 18, color = 'black', hjust = 0.5))





###################################################################################################
###################################################################################################

# RELATEDNESS #

###################################################################################################
###################################################################################################
rel <- read.table("/Users/lanceoconnor/Desktop/JohnsHopkins/classes/AnnotateAgenome/Annotate_a_Genome_Sp2026/exercise_4_GWAS_scripts/RELATE/relatedness.genome", header = TRUE)

# Filter related pairs (PI_HAT > 0.2)
related_pairs <- rel %>%
  filter(PI_HAT > 0.2) %>%
  select(FID1, IID1, FID2, IID2, PI_HAT)

cat("=== RELATED PAIRS (PI_HAT > 0.2) ===\n")
cat("Number of related pairs:", nrow(related_pairs), "\n\n")

if (nrow(related_pairs) > 0) {
  print(related_pairs)
  
  # Strategy: Remove individual with lower call rate or randomly pick one
  # Here we simply remove the second individual from each pair
  # For more sophisticated approach, consider missingness rates
  
  # Get unique individuals to remove (avoid removing both from a pair)
  to_remove <- related_pairs %>%
    select(FID2, IID2) %>%
    distinct() %>%
    rename(FID = FID2, IID = IID2)
  
  # Check if any "to_remove" individuals are also in the "keep" column
  # This avoids removing someone who should be kept
  to_keep <- related_pairs %>%
    select(FID1, IID1) %>%
    distinct() %>%
    rename(FID = FID1, IID = IID1)
  
  # Final removal list: remove IID2 unless they appear as IID1 elsewhere
  final_remove <- to_remove %>%
    anti_join(to_keep, by = c("FID", "IID"))
  
  # If still have pairs, just take one from each pair
  if (nrow(final_remove) == 0 && nrow(related_pairs) > 0) {
    final_remove <- related_pairs %>%
      slice(1) %>%
      select(FID2, IID2) %>%
      rename(FID = FID2, IID = IID2)
  }
  
  cat("\n=== INDIVIDUALS TO REMOVE ===\n")
  print(final_remove)
  
  # Save to file
  write.table(final_remove, "/Users/lanceoconnor/Desktop/JohnsHopkins/classes/AnnotateAgenome/Annotate_a_Genome_Sp2026/exercise_4_GWAS_scripts/RELATE/related_samples_to_remove.txt",
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
} else {
  cat("No related pairs found!\n")
  # Create empty file
  file.create("related_samples_to_remove.txt")
}

# Summary statistics
cat("\n=== PI_HAT SUMMARY ===\n")
cat("Mean PI_HAT:", mean(rel$PI_HAT), "\n")
cat("Max PI_HAT:", max(rel$PI_HAT), "\n")
cat("Pairs with PI_HAT > 0.5:", sum(rel$PI_HAT > 0.5), "\n")
cat("Pairs with PI_HAT > 0.25:", sum(rel$PI_HAT > 0.25), "\n")
cat("Pairs with PI_HAT > 0.2:", sum(rel$PI_HAT > 0.2), "\n")

