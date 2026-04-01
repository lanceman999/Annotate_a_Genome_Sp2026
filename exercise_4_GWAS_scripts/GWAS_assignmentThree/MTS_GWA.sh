#!/bin/bash

plink_files=../PopStrat_assignmentTwo/PLINK_files/QC
outdir=output

# ============================================
# 2a) Multiple testing correction (full dataset)
# ============================================

# --adjust flag adds Bonferroni, FDR, and other corrections
plink --bfile $plink_files/merged_1000G_HapMap \
      --assoc \
      --adjust \
      --out $outdir/assoc_adjusted

# For logistic regression with adjustment
plink --bfile $plink_files/merged_1000G_HapMap \
      --logistic \
      --covar $plink_files/mds_covariates.txt \
      --covar-number 1-10 \
      --adjust \
      --out $outdir/logistic_adjusted

# ============================================
# 2b) Subset SNPs for permutation testing
# ============================================

# Option 1: Random subset of SNPs
# Extract ~10,000 random SNPs for permutation testing
awk '{print $2}' $plink_files/merged_1000G_HapMap.bim | sort -R | head -10000 > $outdir/random_snps.txt

# Create subset PLINK files
plink --bfile $plink_files/merged_1000G_HapMap \
      --extract $outdir/random_snps.txt \
      --make-bed \
      --out $outdir/subset_for_perm

# ============================================
# 2c) Permutation-based association testing
# ============================================

# Run permutation test on subset
# --mperm: number of permutations (e.g., 10000)
plink --bfile $outdir/subset_for_perm \
      --assoc \
      --mperm 10000 \
      --out $outdir/perm_assoc

# ============================================
# Check outputs
# ============================================

echo "=== Adjusted association results ==="
head $outdir/assoc_adjusted.assoc.adjusted
wc -l $outdir/assoc_adjusted.assoc.adjusted

echo ""
echo "=== Permutation results ==="
head $outdir/perm_assoc.assoc.mperm
wc -l $outdir/perm_assoc.assoc.mperm
