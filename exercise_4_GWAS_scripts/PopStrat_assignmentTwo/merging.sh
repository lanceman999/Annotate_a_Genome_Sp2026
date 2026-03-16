#!/bin/bash

outdir=PLINK_files/QC

# ============================================
# 5a) Merge the harmonized datasets
# ============================================

plink --bfile $outdir/1000G_final \
      --bmerge $outdir/HapMap_final \
      --make-bed \
      --out $outdir/merged_1000G_HapMap

# ============================================
# 5b) Report merged dataset counts
# ============================================

echo ""
echo "=== Merged dataset ==="
echo "Variants:"
wc -l $outdir/merged_1000G_HapMap.bim
echo "Samples:"
wc -l $outdir/merged_1000G_HapMap.fam
