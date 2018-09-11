#!/bin/bash

#TISSUE=WholeBlood
#GENE=ENSG00000147852.11
#CHROM=9
#SUMSTATS=/storage/mgymrek/gtex/gwas/summarystats/hgb_N172925_wide_form.tsv.gz
#NAME=VLDLR
#tmpdir=tmp

TISSUE=$1
GENE=$2
CHROM=$3
SUMSTATS=$4
NAME=$5
tmpdir=$(mktemp -d)

# coloc package needs:
# beta
# varbeta

echo $tmpdir

# Extract beta and pvals for gtex/gwas
cat /storage/szfeupe/Runs/650GTEx_estr/Analysis_by_Tissue/${TISSUE}/SNP_Analysis/Lin_Reg_Out | grep -w ${GENE} | cut -f 5,10,11,13 | sort -k1,1g > ${tmpdir}/gtex_locs.txt
start=$(cat ${tmpdir}/gtex_locs.txt | datamash min 1)
end=$(cat ${tmpdir}/gtex_locs.txt | datamash max 1)
zcat ${SUMSTATS} | awk -v "chrom=$CHROM" '($4==chrom)' | \
    awk -v"start=$start" -v "end=$end" '($5>=start && $5<=end)' | \
    cut -f 5,11,12,13,16 | sort -k1,1g > ${tmpdir}/gwas_locs.txt
echo "snppos,gwas.beta,gwas.varbeta,gwas.p,gwas.maf,gtex.beta,gtex.varbeta,gtex.p" | sed 's/,/ /g' > ${NAME}_coloc.txt
join ${tmpdir}/gwas_locs.txt ${tmpdir}/gtex_locs.txt | \
    awk -F" " '{print $1 " " $2 " " $3**2 " " $4 " " $5 " " $6 " " $7**2 " " $8}' >> ${NAME}_coloc.txt
