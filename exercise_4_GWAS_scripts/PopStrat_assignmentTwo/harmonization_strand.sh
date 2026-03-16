#!/bin/bash

outdir=PLINK_files/QC

# ============================================
# 4a) Attempt initial merge to identify problems
# ============================================

# First attempt to merge - this will FAIL but generate a list of problem SNPs
plink --bfile $outdir/1000G_harmonized \
      --bmerge $outdir/HapMap_shared \
      --make-bed \
      --out $outdir/merged_test

# This creates: merged_test-merge.missnp (SNPs with strand issues)

# ============================================
# 4b) Flip strand for mismatched SNPs in 1000G
# ============================================

# If missnp file exists, flip those SNPs and try again
if [ -f $outdir/merged_test-merge.missnp ]; then
    echo "=== Strand mismatches found ==="
    wc -l $outdir/merged_test-merge.missnp
    
    # Flip the problematic SNPs in 1000G
    plink --bfile $outdir/1000G_harmonized \
          --flip $outdir/merged_test-merge.missnp \
          --make-bed \
          --out $outdir/1000G_flipped
    
    # ============================================
    # 4c) Attempt merge again after flipping
    # ============================================
    
    plink --bfile $outdir/1000G_flipped \
          --bmerge $outdir/HapMap_shared \
          --make-bed \
          --out $outdir/merged_test2
fi

# ============================================
# 4d) Remove SNPs that still can't be reconciled
# ============================================

# If there are still problem SNPs after flipping, exclude them
if [ -f $outdir/merged_test2-merge.missnp ]; then
    echo "=== Irreconcilable SNPs ==="
    wc -l $outdir/merged_test2-merge.missnp
    
    # Remove these SNPs from BOTH datasets
    plink --bfile $outdir/1000G_flipped \
          --exclude $outdir/merged_test2-merge.missnp \
          --make-bed \
          --out $outdir/1000G_final
    
    plink --bfile $outdir/HapMap_shared \
          --exclude $outdir/merged_test2-merge.missnp \
          --make-bed \
          --out $outdir/HapMap_final
else
    # No more problem SNPs - just rename
    cp $outdir/1000G_flipped.bed $outdir/1000G_final.bed
    cp $outdir/1000G_flipped.bim $outdir/1000G_final.bim
    cp $outdir/1000G_flipped.fam $outdir/1000G_final.fam
    
    cp $outdir/HapMap_shared.bed $outdir/HapMap_final.bed
    cp $outdir/HapMap_shared.bim $outdir/HapMap_final.bim
    cp $outdir/HapMap_shared.fam $outdir/HapMap_final.fam
fi

# ============================================
# 4e) Report final counts
# ============================================

echo ""
echo "=== Final counts after harmonization ==="
echo "1000G final variants:"
wc -l $outdir/1000G_final.bim
echo "1000G final samples:"
wc -l $outdir/1000G_final.fam

echo ""
echo "HapMap final variants:"
wc -l $outdir/HapMap_final.bim
echo "HapMap final samples:"
wc -l $outdir/HapMap_final.fam
