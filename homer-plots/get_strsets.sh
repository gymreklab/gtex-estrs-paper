#!/bin/bash

# TODO /storage/mgymrek/gtex-estrs/revision/homer-plots/encodedata/peaks/dnase-combined-GM12878-fat-skin-nerve.pm1000.bed
# dnase only look at things close to some dnase peak

OUTDIR=/storage/mgymrek/gtex-estrs/revision/homer-plots/strsets
MASTER=/storage/mgymrek/gtex-estrs/revision/mastertables/
CAUSAL=/storage/mgymrek/gtex-estrs/revision/figures/SuppTable_CAVIAR.tsv
SCORE=0.3 # caviar threshold
DNASE=/storage/mgymrek/gtex-estrs/revision/homer-plots/encodedata/peaks/dnase-combined-GM12878-fat-skin-nerve.pm1000.bed 

# Get all STRs by period
for period in $(seq 1 6)
do
    cat ${MASTER}/*.tab | \
	grep -v chrom | awk -v"period=$period" '(length($9)==period) {print $1 "\t" $3 "\t" $8}'| \
	sort | uniq > ${OUTDIR}/ALLSTRs_period${period}.bed
    intersectBed -a ${OUTDIR}/ALLSTRs_period${period}.bed -b ${DNASE} -u > ${OUTDIR}/ALLSTRs_period${period}_dnase.bed
done

# Get aggregate eSTRs by period
for period in $(seq 1 6)
do
    cat ${CAUSAL} | csvcut -t -c chrom,str.start,str.end,str.motif.forward,score | \
	csvformat -T | grep -v chrom | \
	awk -v"period=$period" -v"score=$SCORE" \
	'(length($4)==period && $5>=score) {print $1 "\t" $2 "\t" $3}' | \
	sort | uniq > ${OUTDIR}/ALLCAUSAL_period${period}.bed
    intersectBed -a ${OUTDIR}/ALLCAUSAL_period${period}.bed -b ${DNASE} -u > ${OUTDIR}/ALLCAUSAL_period${period}_dnase.bed
    cat ${CAUSAL} | csvcut -t -c chrom,str.start,str.end,str.motif.forward,score | \
	csvformat -T | grep -v chrom | \
	awk -v"period=$period" -v"score=$SCORE" \
	'(length($4)==period) {print $1 "\t" $2 "\t" $3}' | \
	sort | uniq > ${OUTDIR}/ALLESTR_period${period}.bed
    intersectBed -a ${OUTDIR}/ALLESTR_period${period}.bed -b ${DNASE} -u > ${OUTDIR}/ALLESTR_period${period}_dnase.bed
done

