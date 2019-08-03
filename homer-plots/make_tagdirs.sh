#!/bin/bash

GENOME=/storage/resources/dbase/human/hg19/hg19.fa
OUTDIR=/storage/mgymrek/gtex-estrs/revision/homer-plots

for f in POLR2A_H1-hESC #POLR2A_GM12878 polr2a-chipseq-heart-ENCFF643EGO polr2a-chipseq-lung-ENCFF719ANC polr2a-chipseq-tibial_nerve-ENCFF750HDH  #dnase-fat-ENCFF880CAD dnase-esophagus-ENCFF208ILP dnase-skin-ENCFF238BRB #DNAse_GM12878 dnase-heart-ENCFF702IJE dnase-lung-ENCFF803ZER dnase-seq-adultheart-ENCFF185ISK dnase-tibial_artery-ENCFF223EKZ dnase-tibial_nerve-ENCFF226ZCG polr2a-chipseq-heart-ENCFF643EGO polr2a-chipseq-lung-ENCFF719ANC polr2a-chipseq-tibial_nerve-ENCFF750HDH
do
    echo makeTagDirectory ${OUTDIR}/tagdirs/${f} ${OUTDIR}/encodedata/${f}.bam \
	-genome ${GENOME}
done | xargs -P10 -n1 -I% sh -c "%"
