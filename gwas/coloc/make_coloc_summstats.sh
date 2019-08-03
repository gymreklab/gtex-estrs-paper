#!/bin/bash

# Make coloc files from all the summary stats files
# Files have chrom,start,end,beta,sebeta,p,maf

OUTDIR=/storage/mgymrek/gtex-estrs/revision/coloc/summstats

# From Glass Lab (http://homer.ucsd.edu/iholtman/Summary_Stats/Summary_Stats.tar.gz)
INT=/storage/mgymrek/gtex-estrs/revision/coloc/summstats/tmp/intelligence/SavageJansen_2018_intelligence_metaanalysis.txt 
cat ${INT} | grep -v UNIQUE | \
    awk '{print $3 "\t" $4 "\t" $4+1 "\t" $9 "\t" $10 "\t" $11 "\t" $7}' | \
    sort -k1,1 -k2,2n | \
    bgzip -c > ${OUTDIR}/Intelligence_coloc.bed.gz
tabix -p bed ${OUTDIR}/Intelligence_coloc.bed.gz

ALZ=/storage/mgymrek/gtex-estrs/revision/coloc/summstats/tmp/Alzheimer/4_UK_Biobank_IGAP_17May2018
cat ${ALZ} | grep -v DIR | \
    awk '{print $8 "\t" $9 "\t" $9+1 "\t" $4 "\t" $5 "\t" $6 "\t" "NA"}' | \
    sort -k1,1 -k2,2n | \
    bgzip -c > ${OUTDIR}/Alzheimers_coloc.bed.gz
tabix -p bed ${OUTDIR}/Alzheimers_coloc.bed.gz

# Height
HT=/storage/mgymrek/gtex/gwas/summarystats//Meta-analysis_Wood_et_al+UKBiobank_2018.txt.gz
zcat ${HT} | grep -v Tested | \
    awk '{print $1 "\t" $2 "\t" $2+1 "\t" $7 "\t" $8 "\t" $9 "\t" $6}' | \
    sort -k1,1 -k2,2n | \
    bgzip -c > ${OUTDIR}/height_coloc.bed.gz
tabix -p bed ${OUTDIR}/height_coloc.bed.gz

# SCZ
SCZ_PGC=/storage/mgymrek/gtex/gwas/summarystats/ckqny.scz2snpres.gz
zcat ${SCZ_PGC} | grep -v hg19chrc | \
    awk '{print $1 "\t" $5 "\t" $5+1 "\t" log($7) "\t" $8/$7 "\t" $9 "\t" "NA"}' | \
    sed 's/chr//' | \
    bgzip -c > ${OUTDIR}/scz_pgc_coloc.bed.gz
tabix -p bed ${OUTDIR}/scz_pgc_coloc.bed.gz

# IBD
IBD=/storage/mgymrek/gtex/gwas/summarystats/EUR.IBD.gwas_info03_filtered.assoc 
cat ${IBD} | grep -v Direction | \
    awk '{print $1 "\t" $3 "\t" $3+1 "\t" log($9) "\t" $10/$9 "\t" $11 "\t" "NA"}' | \
    sort -k1,1 -k2,2n | \
    bgzip -c > ${OUTDIR}/ibd_coloc.bed.gz
tabix -p bed ${OUTDIR}/ibd_coloc.bed.gz
