#!/bin/bash

# Get sample lists for each tissue
TISSUES="Adipose-Subcutaneous Adipose-Visceral Artery-Aorta Artery-Tibial Brain-Caudate Brain-Cerebellum Cells-Transformedfibroblasts Esophagus-Mucosa Esophagus-Muscularis Heart-LeftVentricle Lung Muscle-Skeletal Nerve-Tibial Skin-NotSunExposed Skin-SunExposed Thyroid WholeBlood"
OUTIDR=/storage/mgymrek/gtex-estrs/revision/caviar/

for tissue in $TISSUES
do
    cat /storage/szfeupe/Runs/650GTEx_estr/Analysis_by_Tissue/${tissue}/Corr_Expr.csv  | \
	cut -f 1 -d',' | grep -v ENSG | sed 's/"//g' > ${OUTDIR}/samples/${tissue}.samples
done

# Index genotype data
cp /storage/szfeupe/Runs/650GTEx_estr/Genotypes/NormalizedGenotypes.table /storage/mgymrek/gtex-estrs/revision/genotypes/GTExNormalizedSTRGenotypes.table
bgzip /storage/mgymrek/gtex-estrs/revision/genotypes/GTExNormalizedSTRGenotypes.table
tabix -s 1 -b 2 -e 2 -S 1 /storage/mgymrek/gtex-estrs/revision/genotypes/GTExNormalizedSTRGenotypes.table.gz

for chrom in $(seq 1 22) X
do
    cp /storage/szfeupe/Runs/650GTEx_estr/SNP_Analysis/chr${chrom}.tab /storage/mgymrek/gtex-estrs/revision/genotypes/GTExNormalizedSNPGenotypes_chr${chrom}.table
    bgzip /storage/mgymrek/gtex-estrs/revision/genotypes/GTExNormalizedSNPGenotypes_chr${chrom}.table
    tabix -s 1 -b 2 -e 2 -S 1 /storage/mgymrek/gtex-estrs/revision/genotypes/GTExNormalizedSNPGenotypes_chr${chrom}.table.gz
done

# SNP Genotype data is on aws s3 at s3://gtex-estr/snp_gts_chr${chrom}.tab.gz
# STR Genotype data is on aws s3 at s3://gtex-estr/str_gts.tab.gz
