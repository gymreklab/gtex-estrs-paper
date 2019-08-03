#!/bin/bash

# process ALL genes that we have STR data for

OUTDIR=/storage/mgymrek/gtex-estrs/revision/caviar/batches/
#MASHR=/storage/mgymrek/gtex-estrs/revision/mashr

TISSUES="Adipose-Subcutaneous Adipose-Visceral Artery-Aorta Artery-Tibial Brain-Caudate Brain-Cerebellum Cells-Transformedfibroblasts Esophagus-Mucosa Esophagus-Muscularis Heart-LeftVentricle Lung Muscle-Skeletal Nerve-Tibial Skin-NotSunExposed Skin-SunExposed Thyroid WholeBlood"

BATCHSIZE=100

for tissue in $TISSUES
do
    mkdir -p ${OUTDIR}/${tissue}
    for chrom in $(seq 1 22)
    do
	mkdir -p ${OUTDIR}/${tissue}/chr${chrom}
	cat /storage/mgymrek/gtex-estrs/revision/strreg/${tissue}_strreg.tab  | \
	    grep -w chr${chrom} | \
	    grep -v gene | cut -f 1 | uniq > ${OUTDIR}/${tissue}/chr${chrom}/allgenes.txt
#	cat ${MASHR}/output-strs/sig-bytissue/${tissue}-genelevel-estrs.tsv | \
#	    grep _chr${chrom}_ | cut -f 1 -d'_' | sort | uniq > ${OUTDIR}/${tissue}/chr${chrom}/allgenes.txt
	split -l ${BATCHSIZE} ${OUTDIR}/${tissue}/chr${chrom}/allgenes.txt ${OUTDIR}/${tissue}/chr${chrom}/batch_
    done
done
