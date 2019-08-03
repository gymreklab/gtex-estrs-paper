#!/bin/bash

for chrom in $(seq 1 22)
do
    R CMD BATCH '--args chrom='"${chrom}"'' runMashr_snps.R /storage/mgymrek/gtex-estrs/revision/mashr/output-snps-bychrom/chr${chrom}/mashr.snps_chr${chrom}.log
done
