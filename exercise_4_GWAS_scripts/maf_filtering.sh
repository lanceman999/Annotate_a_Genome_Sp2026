!/bin/bash

workdir="/Users/lanceoconnor/Downloads/lecture_4"
outdir=MAF

plink -bfile $workdir/HapMap_3_r3_1 --freq --out $outdir/allele_freq
