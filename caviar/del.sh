#!/bin/bash

set -e

chrom=3

TISSUES="Adipose-Subcutaneous Adipose-Visceral Artery-Aorta Artery-Tibial Brain-Caudate Brain-Cerebellum Cells-Transformedfibroblasts Esophagus-Mucosa Esophagus-Muscularis Heart-LeftVentricle Lung Muscle-Skeletal Nerve-Tibial Skin-NotSunExposed Skin-SunExposed Thyroid WholeBlood"

for tissue in $TISSUES
do

echo $tissue
# Get old zscores from linreg...
#echo ",${tissue}" | sed 's/,/\t/g' > test.${tissue}.str.old.zscores
#cat /storage/szfeupe/Runs/650GTEx_estr/Analysis_by_Tissue/Artery-Aorta/Lin_Reg_Out | grep ENSG00000163933.5 | awk '{print $1"_"$2"_"$4 "\t" $9/$10}' >> test.${tissue}.str.old.zscores
#echo ",${tissue}" | sed 's/,/\t/g' > test.${tissue}.snp.old.zscores
#cat /storage/szfeupe/Runs/650GTEx_estr/Analysis_by_Tissue/Artery-Aorta/SNP_Analysis/Lin_Reg_Out | grep ENSG00000163933.5 | awk '{print $2"_"$3"_"$5 "\t" $10/$11}' >> test.${tissue}.snp.old.zscores

# Get original zscores from linreg
#echo ",${tissue}" | sed 's/,/\t/g' > test.${tissue}.str.zscores
#cat /storage/mgymrek/gtex-estrs/revision/strreg/${tissue}_strreg.tab | grep -v gene | awk '{print $2"_"$3"_"$5 "\t" ($10/$11)}' | grep ENSG00000163933.5 >> test.${tissue}.str.zscores

#echo ",${tissue}" | sed 's/,/\t/g' > test.${tissue}.snp.zscores
#cat /storage/mgymrek/gtex-estrs/revision/snpreg/${tissue}_snpreg.tab | grep -v gene | awk '{print $1"_"$2"_"$4 "\t" ($6/$7)}' | grep ENSG00000163933.5 >> test.${tissue}.snp.zscores


#    --zstr test.${tissue}.str.zscores \
#    --zsnp test.${tissue}.snp.zscores \

#    --zsnp /storage/mgymrek/gtex-estrs/revision/mashr/output-snps/zscores.tsv \
#    --zstr /storage/mgymrek/gtex-estrs/revision/mashr/output-strs/sig-bytissue/${tissue}-locuslevel-estrs.tsv \

mkdir -p tmp/${tissue}
./RunCaviarGTEx.py \
    --tissue ${tissue} \
    --zstr test.${tissue}.str.zscores --zsnp test.${tissue}.snp.zscores \
    --samples /storage/mgymrek/gtex-estrs/revision/samples/${tissue}.samples \
    --strgt /storage/mgymrek/gtex-estrs/revision/genotypes/GTExNormalizedSTRGenotypes.table.gz \
    --snpgt /storage/mgymrek/gtex-estrs/revision/genotypes/GTExNormalizedSNPGenotypes_chr${chrom}.table.gz \
    --out test.${tissue}.tab \
    --genes ENSG00000163933.5 \
    --tmpdir tmp/${tissue} \
    --num-causal 2 \
    --mingt 3 --recompute-z --use-topn-strs 1 \
    --expr /storage/szfeupe/Runs/650GTEx_estr/Analysis_by_Tissue/Review_Rerun/${tissue}/Corr_Expr.csv
done
