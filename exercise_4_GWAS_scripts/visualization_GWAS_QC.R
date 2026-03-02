library(ggplot2)
library(readr)
library(dplyr)


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

# SNP MISSINGNESS #

###################################################################################################
###################################################################################################
