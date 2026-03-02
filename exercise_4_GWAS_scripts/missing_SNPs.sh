#!/bin/bash

dir="/Users/lanceoconnor/Downloads/lecture_4"

plink --bfile $dir/HapMap_3_r3_1 --missing --out missing_snps/missing_SNPs_out
