#!/bin/bash

set -e

TISSUES="Adipose-Subcutaneous Adipose-Visceral Artery-Aorta Artery-Tibial Brain-Caudate Brain-Cerebellum Cells-Transformedfibroblasts Esophagus-Mucosa Esophagus-Muscularis Heart-LeftVentricle Lung Muscle-Skeletal Nerve-Tibial Skin-NotSunExposed Skin-SunExposed Thyroid WholeBlood"
for tissue in $TISSUES
do
    chrom=3
	echo $tissue $chrom
	grep "ENSG00000163933.5\|ID"  /storage/mgymrek/gtex-estrs/revision/mashr/output-strs/sig-bytissue/${tissue}-genelevel-estrs.tsv > ${tissue}.tmp.tsv
	./RunESTRAnova.py \
	    --sigsnps /storage/mgymrek/gtex-estrs/revision/mashr/output-snps/sig-bytissue/${tissue}-genelevel-estrs.tsv \
	    --sigstrs ${tissue}.tmp.tsv \
	    --samples /storage/mgymrek/gtex-estrs/revision/samples/${tissue}.samples \
	    --strgt /storage/mgymrek/gtex-estrs/revision/genotypes/GTExNormalizedSTRGenotypes.table.gz \
	    --snpgt /storage/mgymrek/gtex-estrs/revision/genotypes/GTExNormalizedSNPGenotypes_chr${chrom}.table.gz \
	    --chrom chr${chrom} \
	    --out ${tissue}.tab \
	    --mingt 3 \
	    --expr /storage/szfeupe/Runs/650GTEx_estr/Analysis_by_Tissue/Review_Rerun/${tissue}/Corr_Expr.csv
done
