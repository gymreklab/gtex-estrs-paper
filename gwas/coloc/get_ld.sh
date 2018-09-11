#!/bin/bash

source params.sh

PREFIX=$1
GENE=$2
TISSUE=$3
STRLOCUS=$4

INFILE=${OUTDIR}/tmp/${PREFIX}_${GENE}_${TISSUE}_combined.txt

head $INFILE
exit 1

SFILE=${OUTDIR}/tmp/${PREFIX}_${GENE}_${TISSUE}_snpstr_loci.txt
OUTFILE=${OUTDIR}/tmp/${PREFIX}_${GENE}_${TISSUE}_LD.txt

/home/mgymrek/workspace/ssc-imputation/snpstr-ld/snp_str_ld_calculator.py \
    --str-vcf ${STRVCF} --snp-vcf ${SNPVCF} --loci-file ${SFILE} \
    --use-info-start --mincount 3 --usefilter --use-gb
