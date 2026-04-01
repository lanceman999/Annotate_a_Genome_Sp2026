#!/bin/bash

# ============================================
# PRS Assignment - Question 4.1
# LD Clumping (Target LD)
# ============================================

# Set directories
workdir=/Users/lanceoconnor/Desktop/JohnsHopkins/classes/AnnotateAgenome/Annotate_a_Genome_Sp2026/lecture4_PRS/output
hapmap_dir=/Users/lanceoconnor/Downloads/lecture_4

cd $workdir

# ============================================
# Step 1: Copy/link HapMap files to working directory
# ============================================

echo "=== Setting up HapMap files ==="

# Create symbolic links to HapMap files
ln -sf $hapmap_dir/HapMap_3_r3_1.bed .
ln -sf $hapmap_dir/HapMap_3_r3_1.bim .
ln -sf $hapmap_dir/HapMap_3_r3_1.fam .

# Verify files exist
ls -la HapMap_3_r3_1.*

# ============================================
# Step 2: Prepare summary stats for PLINK clumping
# ============================================

echo ""
echo "=== Preparing summary stats for clumping ==="

# PLINK --clump expects specific column names: SNP, P
# Our file already has these columns, but let's verify format

head discovery_gwas_sim_sumstats_cleaned.tsv

# ============================================
# Step 3: LD Clumping with PLINK
# ============================================

echo ""
echo "=== Running LD Clumping ==="

# Parameters:
# --clump-p1 1        : Include all SNPs (p-value threshold for index SNPs)
# --clump-p2 1        : Include all SNPs (p-value threshold for clumped SNPs)
# --clump-r2 0.1      : LD r² threshold
# --clump-kb 250      : 250 kb window

plink --bfile HapMap_3_r3_1 \
      --clump discovery_gwas_sim_sumstats_cleaned.tsv \
      --clump-p1 1 \
      --clump-p2 1 \
      --clump-r2 0.1 \
      --clump-kb 250 \
      --out clump

# ============================================
# Step 4: Check output
# ============================================

echo ""
echo "=== Clumping Results ==="

# View clumped output
head clump.clumped
wc -l clump.clumped

# ============================================
# Step 5: Extract clumped SNP list
# ============================================

echo ""
echo "=== Extracting clumped SNPs ==="

# Extract SNP column (column 3) from clumped file
# Skip header, remove empty lines
awk 'NR>1 && $3 != "" {print $3}' clump.clumped > clumped.snps

# Count SNPs
echo "Number of clumped (index) SNPs:"
wc -l clumped.snps

# Preview
echo ""
echo "First 10 clumped SNPs:"
head clumped.snps

echo ""
echo "=== Done! ==="
echo "Output files:"
echo "  1. clump.clumped - Full clumping results"
echo "  2. clumped.snps - List of index SNPs"
