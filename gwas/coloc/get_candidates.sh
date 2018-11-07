#!/bin/bash

# Make output file $OUTDIR/coloc/loci/${prefix}_candidates.tab
# Has columns:
# gene, tissue, rsid, gene.name

source params.sh

TRAIT=$1
PREFIX=$2

GWASHITS=/storage/mgymrek/gtex/gwas/results-Sep2018/str_gwas_ld_COMBINED_eSTR_${TRAIT}.tab
# ENSG00000004838.9Artery-Aortars2073499ZMYND10
cat ${GWASHITS} | awk -v"mincausal=$MINCAUSAL" -F"\t" '($11>mincausal)' | \
    awk -F"\t" '{print $7 "\t" $10 "\t" $13}' | grep -v rsid | \
    ./get_best_tissue.py 3 | ./get_gene.py $ENSGENE 2 | \
    awk '{print $4 "\t" $3 "\t" $2 "\t" $1}' | sort | uniq \
    > $OUTDIR/loci/${PREFIX}_candidates.tab 
#    cut -f 7,9,12 #| ./get_best_tissue.py 3 | ./get_gene.py $ENSGENE 2 | \
#    awk '{print $4 "\t" $3 "\t" $1 "\t" $2}' | sort | uniq > $OUTDIR/loci/${PREFIX}_candidates.tab 

