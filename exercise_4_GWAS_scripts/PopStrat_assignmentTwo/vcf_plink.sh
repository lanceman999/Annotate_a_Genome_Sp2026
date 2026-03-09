#!/bin/bash

workdir="/Users/lanceoconnor/Downloads/lecture_4"
outdir=PLINK_files

mkdir -p $outdir

# 1) Convert VCF to PLINK bed/bim/fam format
plink --vcf $workdir/ALL.2of4intersection.20100804.genotypes.vcf.gz \
      --make-bed \
      --out $outdir/1000G_converted

# 2) Assign IDs to variants with missing variant names (chr:pos format)
plink --bfile $outdir/1000G_converted \
      --set-missing-var-ids @:# \
      --make-bed \
      --out $outdir/1000G_varIDs 
