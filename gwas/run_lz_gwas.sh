#!/bin/bash

# Usage: ./run_lz_gwas.sh gwasfile ldfile snp locus strstart
# Example: ./run_lz_gwas.sh /storage/mgymrek/gtex/gwas/summarystats/hgb_Astle_epacs.tab 9_2622147__WholeBlood_ld.txt rs369552432 rs2219143 VLDLR_HGB strstart

GWASFILE=$1
RSID=$3
RSID2=$4 # use as ref
NAME=$5
WINDOW=100000
LFILE=$2
STRSTART=$6
#LFILE=9_2622147__WholeBlood_ld.txt

line=$(cat ${GWASFILE} | grep -w ${RSID})
chrom=$(echo $line | cut -f 1 -d' ')
pos=$(echo $line | cut -f 2 -d' ')

echo $GWASFILE $RSID
echo ${chrom}:${pos}

head -n 1 pvals.txt > pvals_${NAME}.txt
cat ${GWASFILE} | awk -v"chrom=$chrom" -v"start=$((pos-$WINDOW))" -v"end=$((pos+$WINDOW))" '(($1==chrom) && ($2>start) && ($3<end))' | \
    awk '{print $1 "\t" $2 "\t" $3 "\t" "chr"$1":"$2 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13 "\t" $14}' \
    >> pvals_${NAME}.txt


/storage/resources/source/locuszoom/bin/locuszoom \
    --epacts pvals_${NAME}.txt \
    --ld ${LFILE} \
    --refsnp ${RSID2} \
    --chr ${chrom} --start $((STRSTART-$WINDOW)) --end $((STRSTART+$WINDOW)) \
    --build hg19 \
    --prefix lz/${NAME} --no-cleanup


