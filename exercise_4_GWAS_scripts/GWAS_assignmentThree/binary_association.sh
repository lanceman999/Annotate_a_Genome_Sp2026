#!/bin/bash

plink_files=../PopStrat_assignmentTwo/PLINK_files/QC
outdir=output


# ============================================
# 1a) Basic association test (unadjusted)
# ============================================

# Chi-square test (no covariates)
plink --bfile $plink_files/merged_1000G_HapMap \
      --assoc \
      --out $outdir/basic_assoc

# ============================================
# 1b) Logistic regression adjusted for MDS covariates
# ============================================

plink --bfile $plink_files/merged_1000G_HapMap \
      --logistic \
      --covar $plink_files/mds_covariates.txt \
      --covar-number 1-10 \
      --out $outdir/logistic_assoc

# ============================================
# Check outputs
# ============================================

echo "=== Basic association results ==="
head $outdir/basic_assoc.assoc
wc -l $outdir/basic_assoc.assoc

echo ""
echo "=== Logistic regression results ==="
head $outdir/logistic_assoc.assoc.logistic
wc -l $outdir/logistic_assoc.assoc.logistic
