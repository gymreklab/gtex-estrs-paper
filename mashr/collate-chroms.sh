#!/bin/bash

workdir=$1
prefix=$2

cat ${workdir}/output-${prefix}/chr1/posterior_lfsr*.tsv | head -n 1 > ${workdir}/output-${prefix}/posterior_lfsr.tsv
cat ${workdir}/output-${prefix}/chr1/posterior_betas*.tsv | head -n 1 > ${workdir}/output-${prefix}/posterior_betas.tsv
cat ${workdir}/output-${prefix}/chr1/posterior_beta_ses*.tsv | head -n 1 > ${workdir}/output-${prefix}/posterior_beta_ses.tsv
cat ${workdir}/output-${prefix}/chr1/posterior_log10bf*.tsv | head -n 1 > ${workdir}/output-${prefix}/posterior_log10bf.tsv

for chrom in $(seq 1 22)
do
    cat ${workdir}/output-${prefix}/chr${chrom}/posterior_lfsr*.tsv | grep -v Brain >> ${workdir}/output-${prefix}/posterior_lfsr.tsv
    cat ${workdir}/output-${prefix}/chr${chrom}/posterior_betas*.tsv | grep -v Brain >> ${workdir}/output-${prefix}/posterior_betas.tsv
    cat ${workdir}/output-${prefix}/chr${chrom}/posterior_beta_ses*.tsv | grep -v Brain >> ${workdir}/output-${prefix}/posterior_beta_ses.tsv
    cat ${workdir}/output-${prefix}/chr${chrom}/posterior_log10bf*.tsv | grep -v Brain >> ${workdir}/output-${prefix}/posterior_log10bf.tsv
done
