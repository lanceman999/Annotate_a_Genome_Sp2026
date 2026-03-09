#!/bin/bash

workdir="/Users/lanceoconnor/Downloads/lecture_4"
outdir=HWE

mkdir -p $outdir

plink --bfile $workdir/HapMap_3_r3_1 --hardy --out $outdir/hwe_stats

plink --bfile $workdir/HapMap_3_r3_1 --hwe 1e-6 --make-bed --out $outdir/HapMap_3_r3_1_HWE

plink --bfile $outdir/HapMap_3_r3_1_HWE --hardy --out $outdir/HWE_filtered
