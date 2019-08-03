#!/bin/bash

STRSETS=/storage/mgymrek/gtex-estrs/revision/homer-plots/strsets
NUCL=/storage/mgymrek/gtex-estrs/revision/homer-plots/encodedata/wgEncodeSydhNsomeGm12878Sig.bedGraph
OUTDIR=/storage/mgymrek/gtex-estrs/revision/homer-plots/composite

for period in $(seq 1 6)
do
    /storage/resources/source/homer/bin/annotatePeaks.pl \
	${STRSETS}/ALLCAUSAL_period${period}.bed \
	hg19 -size 1000 -hist 1 -bedGraph ${NUCL} > ${OUTDIR}/ALLCAUSAL_nucl_period${period}.bed
    /storage/resources/source/homer/bin/annotatePeaks.pl \
	${STRSETS}/ALLSTRs_period${period}.bed \
	hg19 -size 1000 -hist 1 -bedGraph ${NUCL} > ${OUTDIR}/ALLSTRs_nucl_period${period}.bed
done
