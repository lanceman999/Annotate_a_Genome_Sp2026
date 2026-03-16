#!/bin/bash

outdir=PLINK_files/QC
workdir="/Users/lanceoconnor/Downloads/lecture_4"

# ============================================
# 6a) Prune SNPs for LD (linkage disequilibrium)
# ============================================

# MDS works better with independent SNPs
# --indep-pairwise: window size, step, r^2 threshold

plink --bfile $outdir/merged_1000G_HapMap \
      --indep-pairwise 50 5 0.2 \
      --out $outdir/pruned_snps

# Extract pruned SNPs
plink --bfile $outdir/merged_1000G_HapMap \
      --extract $outdir/pruned_snps.prune.in \
      --make-bed \
      --out $outdir/merged_pruned

# ============================================
# 6b) Calculate genetic distance matrix (IBS)
# ============================================

plink --bfile $outdir/merged_pruned \
      --genome \
      --out $outdir/merged_ibs

# ============================================
# 6c) Perform MDS analysis
# ============================================

plink --bfile $outdir/merged_pruned \
      --read-genome $outdir/merged_ibs.genome \
      --cluster \
      --mds-plot 10 \
      --out $outdir/merged_mds

# ============================================
# 6d) Report counts
# ============================================

echo "=== MDS analysis ==="
echo "Variants (after LD pruning):"
wc -l $outdir/merged_pruned.bim
echo "Samples:"
wc -l $outdir/merged_pruned.fam
