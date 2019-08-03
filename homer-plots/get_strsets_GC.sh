#!/bin/bash

OUTDIR=/storage/mgymrek/gtex-estrs/revision/homer-plots/strsets
MASTER=/storage/mgymrek/gtex-estrs/revision/mastertables/
CAUSAL=/storage/mgymrek/gtex-estrs/revision/figures/SuppTable_CAVIAR.tsv
SCORE=0.3 # caviar threshold
TSS=/storage/mgymrek/gtex-estrs/revision/misc/promoter_estrs.tab

# Get all STRs
cat ${MASTER}/*.tab | \
    grep -v chrom | awk '{print $1 "\t" $3 "\t" $8}'| \
    intersectBed -a stdin -b $TSS -wa -wb | awk '{print $1 "\t" $2 "\t" $3 "\tx\t0\t" $NF}' | \
    sort | uniq > ${OUTDIR}/ALLSTRs.bed

# Get all FM-eSTRs
cat ${MASTER}/*.tab | \
    grep -v chrom | awk -F"\t" '($26>=0.3) {print $1 "\t" $3 "\t" $8}'| \
    intersectBed -a stdin -b $TSS -wa -wb | awk '{print $1 "\t" $2 "\t" $3 "\tx\t0\t" $NF}' | \
    sort | uniq > ${OUTDIR}/ALL_FMeSTRs.bed

# Get all CCG/CGG STRs and FM-eSTRs
cat ${MASTER}/*.tab | \
    grep -v chrom | grep -w CCG | awk -F"\t" '{print $1 "\t" $3 "\t" $8}'| \
    intersectBed -a stdin -b $TSS -wa -wb | awk '{print $1 "\t" $2 "\t" $3 "\tx\t0\t" $NF}' | \
    sort | uniq > ${OUTDIR}/ALL_CGG_STRs.bed
cat ${MASTER}/*.tab | \
    grep -v chrom | grep -w CCG | awk -F"\t" '($26>=0.3) {print $1 "\t" $3 "\t" $8}'| \
    intersectBed -a stdin -b $TSS -wa -wb | awk '{print $1 "\t" $2 "\t" $3 "\tx\t0\t" $NF}' | \
    sort | uniq > ${OUTDIR}/ALL_CGG_FMeSTRs.bed

# Get all G4 STRs and FM-eSTRs
cat ${MASTER}/*.tab | \
    grep -v chrom | grep -w "CCGGGG\|CGGGCT\|AGCCCG\|CCCG\|ATGCCC\|CGGGG\|CCGGG\|CCTGGG\|ATGGGG\|G\|AGGGGC\|AGGG\|CCCGG\|CCCTG\|CTGGGG\|AGAGGG\|CCCCGG\|CCCCCG\|CCCCG\|ATCCC\|CGGGGG\|CCCT\|AGGGG\|C\|ACCCCC\|AGCCCC" | awk '{print $1 "\t" $3 "\t" $8}'| \
    intersectBed -a stdin -b $TSS -wa -wb | awk '{print $1 "\t" $2 "\t" $3 "\tx\t0\t" $NF}' | \
    sort | uniq > ${OUTDIR}/ALL_G4_STRs.bed
cat ${MASTER}/*.tab | \
    grep -v chrom | grep -w "CCGGGG\|CGGGCT\|AGCCCG\|CCCG\|ATGCCC\|CGGGG\|CCGGG\|CCTGGG\|ATGGGG\|G\|AGGGGC\|AGGG\|CCCGG\|CCCTG\|CTGGGG\|AGAGGG\|CCCCGG\|CCCCCG\|CCCCG\|ATCCC\|CGGGGG\|CCCT\|AGGGG\|C\|ACCCCC\|AGCCCC" | awk -F"\t" '($26>=0.3) {print $1 "\t" $3 "\t" $8}'| \
    intersectBed -a stdin -b $TSS -wa -wb | awk '{print $1 "\t" $2 "\t" $3 "\tx\t0\t" $NF}' | \
    sort | uniq > ${OUTDIR}/ALL_G4_FMeSTRs.bed
