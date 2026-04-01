#!/bin/bash

# ============================================
# PRS Assignment - Question 4.3
# Compute PRS with PLINK
# ============================================

# Set directories
workdir=/Users/lanceoconnor/Desktop/JohnsHopkins/classes/AnnotateAgenome/Annotate_a_Genome_Sp2026/lecture4_PRS/output
covar_file=/Users/lanceoconnor/Desktop/JohnsHopkins/classes/AnnotateAgenome/Annotate_a_Genome_Sp2026/exercise_4_GWAS_scripts/PopStrat_assignmentTwo/PLINK_files/QC/mds_covariates.txt

cd $workdir

# ============================================
# Step 1: Set up required files
# ============================================

echo "=== Setting up files ==="

# Link covariate file
ln -sf $covar_file covar.tsv

# Create phenotype file from fam (if not already exists)
if [ ! -f pheno.tsv ]; then
    echo "Creating pheno.tsv from .fam file..."
    awk 'BEGIN{OFS="\t"; print "FID","IID","PHENO"} {print $1,$2,$6}' \
        HapMap_3_r3_1.fam > pheno.tsv
fi

# Check files exist
echo ""
echo "Checking required files:"
ls -la HapMap_3_r3_1.bed
ls -la HapMap_3_r3_1.bim
ls -la HapMap_3_r3_1.fam
ls -la pheno.tsv
ls -la covar.tsv
ls -la prs_weights_harmonized_qc.tsv

echo ""
echo "SNP threshold files:"
ls -la prs_*.snps

# ============================================
# Step 2: Prepare score file for PLINK
# ============================================

echo ""
echo "=== Preparing score file ==="

# PLINK --score expects: SNP, A1, BETA
# Our prs_weights_harmonized_qc.tsv has: SNP, A1, BETA, P
# Keep header for --score header option

cut -f1,2,3 prs_weights_harmonized_qc.tsv > prs_score_file.txt

echo "Score file created:"
head prs_score_file.txt
echo "..."
wc -l prs_score_file.txt

# ============================================
# Step 3: Compute PRS for each threshold
# ============================================

echo ""
echo "=== Computing PRS for each threshold ==="

# Define thresholds (must match filenames)
thresholds=("5e-8" "1e-6" "1e-4" "1e-3" "1e-2" "0.05" "0.1" "0.5" "1")

for thresh in "${thresholds[@]}"; do
    
    echo ""
    echo "--- Threshold: P ≤ $thresh ---"
    
    snp_file="prs_${thresh}.snps"
    
    # Check SNP file exists and has content
    if [ ! -s "$snp_file" ]; then
        echo "Warning: $snp_file is empty or missing, skipping..."
        continue
    fi
    
    n_snps=$(wc -l < "$snp_file")
    echo "SNPs in threshold: $n_snps"
    
    # Compute PRS using PLINK
    plink --bfile HapMap_3_r3_1 \
          --score prs_score_file.txt 1 2 3 header \
          --extract "$snp_file" \
          --out "prs_${thresh}"
    
    # Check output
    if [ -f "prs_${thresh}.profile" ]; then
        echo "Created: prs_${thresh}.profile"
        head -3 "prs_${thresh}.profile"
    else
        echo "Error: prs_${thresh}.profile not created!"
    fi
    
done

# ============================================
# Step 4: Summary
# ============================================

echo ""
echo "============================================"
echo "         PRS COMPUTATION SUMMARY           "
echo "============================================"

printf "%-20s %10s %10s\n" "Threshold" "SNPs" "Individuals"
printf "%-20s %10s %10s\n" "---------" "----" "-----------"

for thresh in "${thresholds[@]}"; do
    snp_file="prs_${thresh}.snps"
    profile="prs_${thresh}.profile"
    
    if [ -f "$snp_file" ] && [ -f "$profile" ]; then
        n_snps=$(wc -l < "$snp_file")
        n_ind=$(tail -n +2 "$profile" | wc -l)
        printf "%-20s %10d %10d\n" "P ≤ $thresh" "$n_snps" "$n_ind"
    fi
done

echo "============================================"
echo ""
echo "Profile file columns:"
echo "  FID    - Family ID"
echo "  IID    - Individual ID"
echo "  PHENO  - Phenotype"
echo "  CNT    - Number of non-missing SNPs"
echo "  CNT2   - Number of named alleles"
echo "  SCORE  - PRS (sum of BETA * dosage)"
echo ""
echo "=== Done! ==="
