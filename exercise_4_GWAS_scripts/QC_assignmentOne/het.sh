#!/bin/bash

workdir="/Users/lanceoconnor/Downloads/lecture_4"
outdir=heterozygosity

mkdir -p $outdir

#plink --bfile $workdir/HapMap_3_r3_1 --het --out $outdir/het_stats

#plink --bfile $workdir/HapMap_3_r3_1 --remove $outdir/het_outliers.txt --make-bed --out $outdir/HapMap_3_r3_1_HETfilt

plink --bfile $outdir/HapMap_3_r3_1_HETfilt --het --out $outdir/het_stats_filtered
