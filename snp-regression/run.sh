#!/bin/bash

TISSUE=$1
CHROM=$2

# Get all the files we need from S3
aws s3 cp s3://gtex-estr/snp_gts_chr${CHROM}.tab.gz /scratch/snp_gts_chr${CHROM}.tab.gz
aws s3 cp s3://gtex-estr/snp_gts_chr${CHROM}.tab.gz.tbi /scratch/snp_gts_chr${CHROM}.tab.gz.tbi
aws s3 cp s3://gtex-estr/gencode_gene_annotations_hg19.csv /scratch/gencode_gene_annotations_hg19.csv
aws s3 cp s3://gtex-estr/Corr_Expr_${TISSUE}.csv /scratch/Corr_Expr_${TISSUE}.csv

# Run regression analysis
./LinRegAssociationTest_v2.py \
    --chrom ${CHROM} \
    --expr /scratch/Corr_Expr_${TISSUE}.csv \
    --exprannot /scratch/gencode_gene_annotations_hg19.csv \
    --strgt /scratch/snp_gts_chr${CHROM}.tab.gz \
    --distfromgene 100000 \
    --norm \
    --out Lin_Reg_Out_SNPs_${TISSUE}_${CHROM}.tab \
    --tmpdir /tmp

# Upload results
aws s3 cp Lin_Reg_Out_SNPs_${TISSUE}_${CHROM}.tab s3://gtex-estr/RegressionResults/Lin_Reg_Out_SNPs_${TISSUE}_${CHROM}.tab
