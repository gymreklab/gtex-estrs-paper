#!/bin/bash

source params.sh

PREFIX=$1
GENE=$2
TISSUE=$3
STRLOCUS=$4

INFILE=${OUTDIR}/tmp/${PREFIX}_${GENE}_${TISSUE}_combined.txt
SFILE=${OUTDIR}/tmp/${PREFIX}_${GENE}_${TISSUE}_snpstr_loci.txt
OUTFILE=${OUTDIR}/tmp/${PREFIX}_${GENE}_${TISSUE}_LD.txt

cat $INFILE | grep -v snppos | awk -F" " -v"strloc=$STRLOCUS" '{print strloc "\t" $1}' | sed 's/:/\t/' > ${SFILE}

echo "snppos,ld" | sed 's/,/\t/g' > $OUTFILE
/home/mgymrek/workspace/ssc-imputation/snpstr-ld/snp_str_ld_calculator.py \
    --str-vcf ${STRVCF} --snp-vcf ${SNPVCF} --loci-file ${SFILE} \
    --use-info-start --mincount 3 --usefilter --use-gb | \
    cut -f 2,7 | cut -f 2- -d':' | grep -v locus >> $OUTFILE
