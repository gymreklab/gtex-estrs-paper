#!/bin/bash

set -e

# Clean up old anova runs
rm /storage/mgymrek/gtex-estrs/revision/anova/*.tab

TISSUES="Adipose-Subcutaneous Adipose-Visceral Artery-Aorta Artery-Tibial Brain-Caudate Brain-Cerebellum Cells-Transformedfibroblasts Esophagus-Mucosa Esophagus-Muscularis Heart-LeftVentricle Lung Muscle-Skeletal Nerve-Tibial Skin-NotSunExposed Skin-SunExposed Thyroid WholeBlood"
for tissue in $TISSUES
do
    ZSNP=/storage/mgymrek/gtex-estrs/revision/anova/top-zscores/${tissue}.topz.snps.tsv
    ZSTR=/storage/mgymrek/gtex-estrs/revision/anova/top-zscores/${tissue}.topz.strs.tsv
    # Get top Z scores for SNPs/STRs. Need ID (gene_chr_start) and tissue as columns
    echo "ID,${tissue}" | sed 's/,/\t/g' > ${ZSNP}
    cat /storage/mgymrek/gtex-estrs/revision/snpreg/${tissue}_snpreg.tab | grep -v gene | \
	awk '{print $1 "\t" $1"_"$2"_"$4 "\t" ($6/$7) "\t" (($6>0)?($6/$7):(-1*$6/$7)) }' | \
	datamash -g 1 max 4 -f | cut -f 2,3 >> ${ZSNP}
    echo "ID,${tissue}" | sed 's/,/\t/g' > ${ZSTR}
    cat /storage/mgymrek/gtex-estrs/revision/strreg/${tissue}_strreg.tab | grep -v gene | \
	awk '{print $1 "\t" $1"_"$2"_"$4 "\t" ($6/$7) "\t" (($6>0)?($6/$7):(-1*$6/$7)) }' | \
	datamash -g 1 max 4 -f | cut -f 2,3 >> ${ZSTR}
    for chrom in $(seq 1 22)
    do
	echo ./RunESTRAnova.py \
	    --sigsnps ${ZSNP} \
	    --sigstrs ${ZSTR} \
	    --samples /storage/mgymrek/gtex-estrs/revision/samples/${tissue}.samples \
	    --strgt /storage/mgymrek/gtex-estrs/revision/genotypes/GTExNormalizedSTRGenotypes.table.gz \
	    --snpgt /storage/mgymrek/gtex-estrs/revision/genotypes/GTExNormalizedSNPGenotypes_chr${chrom}.table.gz \
	    --chrom chr${chrom} \
	    --mingt 3 \
	    --out /storage/mgymrek/gtex-estrs/revision/anova/${tissue}_chr${chrom}_anova.tab \
	    --expr /storage/szfeupe/Runs/650GTEx_estr/Analysis_by_Tissue/Review_Rerun/${tissue}/Corr_Expr.csv
    done
done | xargs -n1 -I% -P20 sh -c "%"
