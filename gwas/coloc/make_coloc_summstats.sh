#!/bin/bash

source params.sh
# Make coloc files from all the summary stats files
# Files have chrom,start,end,beta,sebeta,p,maf

# SCZ
SCZ_PGC=/storage/mgymrek/gtex/gwas/summarystats/ckqny.scz2snpres.gz
zcat ${SCZ_PGC} | grep -v hg19chrc | \
    awk '{print $1 "\t" $5 "\t" $5+1 "\t" log($7) "\t" $8/$7 "\t" $9 "\t" "NA"}' | \
    sed 's/chr//' | \
    bgzip -c > ${OUTDIR}/scz_pgc_coloc.bed.gz
tabix -p bed ${OUTDIR}/scz_pgc_coloc.bed.gz

exit 1

# Blood traits
HGB_ASTLE=/storage/mgymrek/gtex/gwas/summarystats/hgb_N172925_wide_form.tsv.gz
BASOP_ASTLE=/storage/mgymrek/gtex/gwas/summarystats/baso_p_gran_N170223_wide_form.tsv.gz
RBC_ASTLE=/storage/mgymrek/gtex/gwas/summarystats/rbc_N172952_wide_form.tsv.gz

for btf in ${BASOP_ASTLE} ${RBC_ASTLE} ${HGB_ASTLE}
do
    prefix=$(basename ${btf} | sed 's/_wide_form.tsv.gz//')
    echo $prefix
    zcat ${btf} | awk '{print $4 "\t" $5 "\t" $5+1 "\t" $11 "\t" $12 "\t" $13 "\t" $16}' | grep -v  EFFECT |\
	bgzip -c > ${OUTDIR}/${prefix}_coloc.bed.gz
    tabix -p bed ${OUTDIR}/${prefix}_coloc.bed.gz
done
