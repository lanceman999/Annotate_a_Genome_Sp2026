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

ggplot(maf, aes(x = MAF)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 0.05, color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "Minor Allele Frequency Distribution",
       x = "MAF",
       y = "Count",
       caption = "Red line = MAF threshold (0.05)") +
  theme_minimal()





