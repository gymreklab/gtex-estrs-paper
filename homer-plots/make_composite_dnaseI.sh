#!/bin/bash

# TODO not done
STRSETS=/storage/mgymrek/gtex-estrs/revision/homer-plots/strsets
OUTDIR=/storage/mgymrek/gtex-estrs/revision/homer-plots/composite
TAGDIR=/storage/mgymrek/gtex-estrs/revision/homer-plots/tagdirs
DNASE="${TAGDIR}/dnase-fat-ENCFF880CAD ${TAGDIR}/dnase-esophagus-ENCFF208ILP ${TAGDIR}/dnase-skin-ENCFF238BRB ${TAGDIR}/DNAse_GM12878 ${TAGDIR}/dnase-heart-ENCFF702IJE ${TAGDIR}/dnase-lung-ENCFF803ZER ${TAGDIR}/dnase-seq-adultheart-ENCFF185ISK ${TAGDIR}/dnase-tibial_artery-ENCFF223EKZ ${TAGDIR}/dnase-tibial_nerve-ENCFF226ZCG"

TISSUES="Adipose-Subcutaneous Adipose-Visceral Artery-Aorta Artery-Tibial Brain-Caudate Brain-Cerebellum Cells-Transformedfibroblasts Esophagus-Mucosa Esophagus-Muscularis Heart-LeftVentricle Lung Muscle-Skeletal Nerve-Tibial Skin-NotSunExposed Skin-SunExposed Thyroid WholeBlood"

for period in $(seq 1 6)
do
    cmd1="/storage/resources/source/homer/bin/annotatePeaks.pl \
	${STRSETS}/ALLCAUSAL_period${period}.bed \
	hg19 -size 1000 -hist 1 -d ${DNASE} > ${OUTDIR}/ALLCAUSAL_dnase_period${period}_nopeak.bed"
    cmd2="/storage/resources/source/homer/bin/annotatePeaks.pl \
	${STRSETS}/ALLSTRs_period${period}.bed \
	hg19 -size 1000 -hist 1 -d ${DNASE} > ${OUTDIR}/ALLSTRs_dnase_period${period}_nopeak.bed"
    echo $cmd1
    echo $cmd2
#    for tissue in $TISSUES
#    do
#	cmd3="/storage/resources/source/homer/bin/annotatePeaks.pl \
#	    ${STRSETS}/${tissue}_period${period}.bed \
#	    hg19 -size 1000 -ghist -hist 1 -d ${DNASE} > ${OUTDIR}/${tissue}_dnase_period${period}.bed"
#	echo $cmd3
#   done
done | xargs -P20 -n1 -I% sh -c "%"
