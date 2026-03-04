!/bin/bash

workdir="/Users/lanceoconnor/Downloads/lecture_4"
outdir=sex_discrep

mkdir -p $outdir

echo "Step 1: Checking sex discrepancies"
plink --bed $workdir/HapMap_3_r3_1.bed --bim $workdir/HapMap_3_r3_1.bim --fam $workdir/HapMap_3_r3_1.fam \
  --check-sex --out $outdir/sex_check

echo "Step 2: Extracting flagged samples"
awk '$5 == "PROBLEM" {print $1, $2}' $outdir/sex_check.sexcheck > $outdir/problem_samples.txt
echo "Found $(wc -l < problem_samples.txt) discrepant samples:"
cat $outdir/problem_samples.txt

echo "Step 3: Imputing sex"
awk 'NR > 1 {
  if ($6 < 0.2) sex = 2;
  else if ($6 > 0.8) sex = 1;
  else sex = 0;
  print $1, $2, sex
}' $outdir/sex_check.sexcheck > $outdir/imputed_sex.txt

plink --bed $workdir/HapMap_3_r3_1.bed --bim $workdir/HapMap_3_r3_1.bim --fam $workdir/HapMap_3_r3_1.fam \
  --update-sex $outdir/imputed_sex.txt \
  --make-bed --out $outdir/HapMap_3_r3_1_sex_corrected
