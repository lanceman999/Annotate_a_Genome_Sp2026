#!/bin/bash

#SBATCH -J filtVCF                   # Job name
#SBATCH -A eande106                     # Allocation name
#SBATCH -p parallel                     # Partition/Queue name
#SBATCH -t 48:00:00                     # Job walltime/duration (hh:mm:ss)
#SBATCH -N 1                            # Number of nodes
#SBATCH -c 24                           # Number of cores

# VCF with 684 isotype reference strains
#input="/vast/eande106/data/c_elegans/WI/variation/20250625/vcf/WI.20250625.hard-filter.isotype.vcf.gz"

#bcftools view -m2 -M2 -e 'CHROM=="MtDNA"' -v snps -q 0.05 -O z -o WI.20250625.hard-filter.isotype.vcf.biallelicSNPs.MAF5.gz $input

# VCF with all 1955 strains
input="/vast/eande106/data/c_elegans/WI/variation/20250625/vcf/WI.20250625.hard-filter.vcf.gz"

bcftools view -m2 -M2 -e 'CHROM=="MtDNA"' -v snps -q 0.05 -O z -o WI.20250625.hard-filter.allStrains.biallelicSNPs.MAF5.gz $input
