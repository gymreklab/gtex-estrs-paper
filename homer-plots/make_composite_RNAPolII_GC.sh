#!/bin/bash

STRSETS=/storage/mgymrek/gtex-estrs/revision/homer-plots/strsets
NUCL=/storage/mgymrek/gtex-estrs/revision/homer-plots/encodedata/wgEncodeSydhNsomeGm12878Sig.bedGraph
OUTDIR=/storage/mgymrek/gtex-estrs/revision/homer-plots/composite
TAGDIR=/storage/mgymrek/gtex-estrs/revision/homer-plots/tagdirs

POL2="${TAGDIR}/POLR2A_GM12878 ${TAGDIR}/polr2a-chipseq-heart-ENCFF643EGO ${TAGDIR}/polr2a-chipseq-lung-ENCFF719ANC ${TAGDIR}/polr2a-chipseq-tibial_nerve-ENCFF750HDH ${TAGDIR}/POLR2A_H1-hESC"

cd /storage/mgymrek/del

for prefix in ALLSTRs ALL_FMeSTRs ALL_CGG_STRs ALL_CGG_FMeSTRs ALL_G4_STRs ALL_G4_FMeSTRs
do
    nohup sh -c "/storage/resources/source/homer/bin/annotatePeaks.pl \
	${STRSETS}/${prefix}.bed \
	hg19 -size 10000 -hist 5 -d ${POL2} > ${OUTDIR}/${prefix}_pol2.bed" &
done
