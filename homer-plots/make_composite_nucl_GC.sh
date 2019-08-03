#!/bin/bash

STRSETS=/storage/mgymrek/gtex-estrs/revision/homer-plots/strsets
NUCL=/storage/mgymrek/gtex-estrs/revision/homer-plots/encodedata/wgEncodeSydhNsomeGm12878Sig.bedGraph
OUTDIR=/storage/mgymrek/gtex-estrs/revision/homer-plots/composite

cd /storage/mgymrek/del

for prefix in ALLSTRs ALL_FMeSTRs ALL_CGG_STRs ALL_CGG_FMeSTRs ALL_G4_STRs ALL_G4_FMeSTRs
do
    nohup sh -c "/storage/resources/source/homer/bin/annotatePeaks.pl \
	${STRSETS}/${prefix}.bed \
	hg19 -size 10000 -hist 5 -bedGraph ${NUCL} > ${OUTDIR}/${prefix}_nucl.bed" &
done
