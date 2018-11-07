#!/bin/bash

source params.sh
PREFIX=$1
GPREFIX=$2

# Get file with:
# p(H0) through p(H4) from coloc
# causality score of the STR
GWASHITS=/storage/mgymrek/gtex/gwas/results-Sep2018/str_gwas_ld_COMBINED_eSTR_${GPREFIX}.tab

echo "chrom,gene,rsid,str.start,tissue,score,ld,p0,p1,p2,p3,p4" | sed 's/,/\t/g'

for resfile in $(ls ${OUTDIR}/${PREFIX}_*_results.txt)
do
    gene=$(basename $resfile | awk -F"_" '{print $(NF-1)}')
    line=$(cat $resfile | head -n 2 | tail -n 1)
    p0=$(echo $line | awk '{print $1}')
    p1=$(echo $line | awk '{print $2}')
    p2=$(echo $line | awk '{print $3}')
    p3=$(echo $line | awk '{print $4}')
    p4=$(echo $line | awk '{print $5}')

    # Get tissue
    tissue=$(cat $OUTDIR/loci/${PREFIX}_candidates.tab | grep -w $gene | cut -f 2 | head -n 1)

    # Get SNP, LD, and STR
    rsid=$(cat $OUTDIR/loci/${PREFIX}_candidates.tab | grep -w $gene | cut -f 4 | head -n 1)
    gwasline=$(cat $GWASHITS | grep -w $rsid | grep -w $gene | sort -k11,11g | tail -n 1)
    str=$(echo $gwasline | cut -f 2 -d' ')
    ld=$(echo $gwasline | cut -f 3 -d' ')
    chrom=$(echo $gwasline | cut -f 1 -d' ')

    # Get score
    score=$(cat $CAUSAL | grep $gene | cut -d',' -f 9 | cut -f 1 -d' ' | head -n 1 )
    echo $chrom,$gene,$rsid,$str,$tissue,$score,$ld,$p0,$p1,$p2,$p3,$p4 | sed 's/,/\t/g'
done
