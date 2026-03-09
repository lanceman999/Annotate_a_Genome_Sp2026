#!/bin/bash

workdir="/Users/lanceoconnor/Downloads/lecture_4"
outdir=RELATE

mkdir -p $outdir


# LD pruning: window 50 SNPs, step 5, r² threshold 0.2
#plink --bfile $workdir/HapMap_3_r3_1 --indep-pairwise 50 5 0.2 --out $outdir/ld_pruned

# Extract pruned SNPs
#plink --bfile $workdir/HapMap_3_r3_1 --extract $outdir/ld_pruned.prune.in --make-bed --out $outdir/HapMap_3_r3_1_pruned

# Compute relatedness
#plink --bfile $outdir/HapMap_3_r3_1_pruned --genome --out $outdir/relatedness

plink --bfile $workdir/HapMap_3_r3_1 --remove $outdir/related_samples_to_remove.txt --make-bed --out $outdir/HapMap_3_r3_1_UNRELATED
