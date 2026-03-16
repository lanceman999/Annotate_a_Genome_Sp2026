#!/bin/bash

outdir=PLINK_files/QC
workdir="/Users/lanceoconnor/Downloads/lecture_4"

# ============================================
# 3a) Create position update file from HapMap
# ============================================

awk '{print $2, $4}' $outdir/HapMap_shared.bim > $outdir/hapmap_positions.txt

# ============================================
# 3b) Update 1000G positions to match HapMap
# ============================================

plink --bfile $outdir/1000G_shared \
      --update-map $outdir/hapmap_positions.txt \
      --make-bed \
      --out $outdir/1000G_harmonized

# ============================================
# 3c) Report and verify
# ============================================

echo "=== After harmonization ==="
echo "Variants:"
wc -l $outdir/1000G_harmonized.bim
echo "Samples:"
wc -l $outdir/1000G_harmonized.fam

echo ""
echo "=== Verify positions now match ==="
echo "1000G harmonized:"
head -5 $outdir/1000G_harmonized.bim
echo ""
echo "HapMap:"
head -5 $outdir/HapMap_shared.bim
