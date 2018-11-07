#!/bin/bash

source params.sh

# chrom gene rsid str.start tissue score ld p0 p1 p2 p3 p4
PREFIX=$1
SUMMFILE=${OUTDIR}/${PREFIX}_coloc_summary.tab
cat $SUMMFILE | grep -v gene | awk '($12>0.5)' > $SUMMFILE.pass

while IFS='' read -r line || [[ -n "$line" ]]; do
    tissue=$(echo $line | cut -f 5 -d' ')
    gene=$(echo $line | cut -f 2 -d' ')
    chrom=$(echo $line | cut -f 1 -d' ')
    start=$(echo $line | cut -f 4 -d' ')
    if [ x"$tissue" == x"tissue" ]; then continue; fi
    ./plot_expr.py $tissue $gene $chrom $start exprplots $PREFIX
done < "$SUMMFILE.pass"
