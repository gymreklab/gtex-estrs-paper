#!/bin/bash
set -e

# See setup.sh for getting files ready

# Clean up
rm -rf /storage/mgymrek/gtex-estrs/revision/caviar/batches/*
rm -rf /storage/mgymrek/gtex-estrs/revision/caviar/tmp/*
rm -f /storage/mgymrek/gtex-estrs/revision/caviar/output/*

# See get_batches.sh for getting gene batches
./get_batches.sh

# Run on batches
TISSUES="Adipose-Subcutaneous Adipose-Visceral Artery-Aorta Artery-Tibial Brain-Caudate Brain-Cerebellum Cells-Transformedfibroblasts Esophagus-Mucosa Esophagus-Muscularis Heart-LeftVentricle Lung Muscle-Skeletal Nerve-Tibial Skin-NotSunExposed Skin-SunExposed Thyroid WholeBlood"

for tissue in $TISSUES
do
    # Set up zscores - use all SNPs
    ZSNP=/storage/mgymrek/gtex-estrs/revision/caviar/zscores/${tissue}.zscores.snps.tsv
    echo ",${tissue}" | sed 's/,/\t/g' > ${ZSNP}
    cat /storage/mgymrek/gtex-estrs/revision/snpreg/${tissue}_snpreg.tab | grep -v gene | \
	awk '{print $1"_"$2"_"$4 "\t" ($6/$7)}' >> ${ZSNP}

    # Set up zscores - use only top STRs per gene, based on original z scores
    ZSTR=/storage/mgymrek/gtex-estrs/revision/caviar/zscores/${tissue}.zscores.strs.tsv
    echo ",${tissue}" | sed 's/,/\t/g' > ${ZSTR}
    cat /storage/mgymrek/gtex-estrs/revision/strreg/${tissue}_strreg.tab | grep -v gene | \
	awk '{print $1 "\t" $1"_"$2"_"$4 "\t" ($6/$7) "\t" (($6>0)?($6/$7):(-1*$6/$7)) }' | \
	datamash -g 1 max 4 -f | cut -f 2,3 >> ${ZSTR}

    # Run CAVIAR per chrom. Use top STR + all SNPs
    for chrom in $(seq 1 22)
    do
	batches=$(ls /storage/mgymrek/gtex-estrs/revision/caviar/batches/${tissue}/chr${chrom}/batch* | awk -F"/" '{print $NF}')
	for batch in $batches
	do
	    mkdir -p /storage/mgymrek/gtex-estrs/revision/caviar/tmp/${tissue}/${chrom}-${batch}/
	    echo ./RunCaviarGTEx.py \
		--zsnp ${ZSNP} --zstr ${ZSTR} \
		--tissue ${tissue} \
		--samples /storage/mgymrek/gtex-estrs/revision/samples/${tissue}.samples \
		--strgt /storage/mgymrek/gtex-estrs/revision/genotypes/GTExNormalizedSTRGenotypes.table.gz \
		--snpgt /storage/mgymrek/gtex-estrs/revision/genotypes/GTExNormalizedSNPGenotypes_chr${chrom}.table.gz \
		--out /storage/mgymrek/gtex-estrs/revision/caviar/output/${tissue}_${chrom}_${batch}.tab \
		--genes-file /storage/mgymrek/gtex-estrs/revision/caviar/batches/${tissue}/chr${chrom}/${batch} \
		--tmpdir /storage/mgymrek/gtex-estrs/revision/caviar/tmp/${tissue}/${chrom}-${batch}/ \
		--num-causal 2 --mingt 3 --recompute-z --use-topn-strs 1 \
		--zthresh 1.96 \
		--expr /storage/szfeupe/Runs/650GTEx_estr/Analysis_by_Tissue/Review_Rerun/${tissue}/Corr_Expr.csv
	done
    done
done | xargs -n1 -I% -P20 sh -c "%"
