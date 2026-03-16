#!/bin/bash

workdir="/Users/lanceoconnor/Downloads/lecture_4"
files=PLINK_files
outdir=PLINK_files/QC

# ============================================
# 2a) QC filtering on 1000 Genomes data
# ============================================

# Filter for:
# --maf 0.01        : Minor allele frequency > 1%
# --geno 0.10       : SNP missingness < 10% (exclude SNPs missing in >10% of samples)
# --hwe 1e-6        : Hardy-Weinberg equilibrium p-value threshold (optional but recommended)

plink --bfile $files/1000G_varIDs \
      --maf 0.01 \
      --geno 0.1 \
      --hwe 1e-6 \
      --make-bed \
      --out $outdir/1000G_QC

# ============================================
# 2b) Report variants and samples after QC
# ============================================

echo "=== After QC ==="
echo "Variants:"
wc -l $outdir/1000G_QC.bim
echo "Samples:"
wc -l $outdir/1000G_QC.fam

# ============================================
# 2c) Find shared SNPs between 1000G and HapMap
# ============================================

# Extract SNP IDs from both datasets
cut -f2 $outdir/1000G_QC.bim > $outdir/1000G_snps.txt
cut -f2 $workdir/HapMap_3_r3_1.bim > $outdir/HapMap_snps.txt

# Find overlapping SNPs
comm -12 <(sort $outdir/1000G_snps.txt) <(sort $outdir/HapMap_snps.txt) > $outdir/shared_snps.txt

echo "=== Shared SNPs ==="
wc -l $outdir/shared_snps.txt

# ============================================
# 2d) Extract shared SNPs from both datasets
# ============================================

# i) Filter 1000G to shared SNPs only
plink --bfile $outdir/1000G_QC \
      --extract $outdir/shared_snps.txt \
      --make-bed \
      --out $outdir/1000G_shared

# ii) Filter HapMap to shared SNPs only
plink --bfile $workdir/HapMap_3_r3_1 \
      --extract $outdir/shared_snps.txt \
      --make-bed \
      --out $outdir/HapMap_shared

# ============================================
# Report final counts
# ============================================

echo ""
echo "=== Final counts ==="
echo "1000G shared variants:"
wc -l $outdir/1000G_shared.bim
echo "1000G shared samples:"
wc -l $outdir/1000G_shared.fam

echo ""
echo "HapMap shared variants:"
wc -l $outdir/HapMap_shared.bim
echo "HapMap shared samples:"
wc -l $outdir/HapMap_shared.fam
