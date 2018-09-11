#!/bin/bash

# Make output file $OUTDIR/coloc/loci/${prefix}_candidates.tab
# Has columns:
# gene, tissue, rsid, gene.name

source params.sh

TRAIT=$1
PREFIX="$2"

cat ${GWASHITS} | grep -w "${TRAIT}" | awk -v"mincausal=$MINCAUSAL" -F"\t" '($10>mincausal)' | \
    cut -f 7,9,12 | ./get_best_tissue.py 3 | ./get_gene.py $ENSGENE 2 | \
    awk '{print $4 "\t" $3 "\t" $1 "\t" $2}' | sort | uniq > $OUTDIR/loci/${PREFIX}_candidates.tab 

