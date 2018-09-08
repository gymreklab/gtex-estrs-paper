#!/bin/bash

# Usage: ./run_lz_gwas.sh gwasfile ldfile snp locus
# Example: ./run_lz_gwas.sh /storage/mgymrek/gtex/gwas/summarystats/hgb_Astle_epacs.tab 9_2622147__WholeBlood_ld.txt rs369552432 rs2219143 VLDLR_HGB

GWASFILE=$1
RSID=$3
RSID2=$4 # use as ref is top is indel
NAME=$5
WINDOW=100000
LFILE=$2
#LFILE=9_2622147__WholeBlood_ld.txt

line=$(cat ${GWASFILE} | grep ${RSID})
chrom=$(echo $line | cut -f 1 -d' ')
pos=$(echo $line | cut -f 2 -d' ')

echo $GWASFILE $RSID
echo ${chrom}:${pos}

head -n 1 pvals.txt > pvals_${NAME}.txt
cat ${GWASFILE} | awk -v"chrom=$chrom" -v"start=$((pos-$WINDOW))" -v"end=$((pos+$WINDOW))" '(($1==chrom) && ($2>start) && ($3<end))' >> pvals_${NAME}.txt

/storage/resources/source/locuszoom/bin/locuszoom \
    --epacts pvals_${NAME}.txt \
    --ld ${LFILE} \
    --chr ${chrom} --start $((pos-$WINDOW)) --end $((pos+$WINDOW)) \
    --build hg19  --source 1000G_March2012 --pop EUR \
    --prefix lz/${NAME} --no-cleanup


