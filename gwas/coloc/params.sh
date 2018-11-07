#!/bin/bash

OUTDIR=/storage/mgymrek/gtex/gwas/summarystats/coloc
#GWASHITS=/storage/mgymrek/gtex/gwas/results-Sep2018/str_gwas_ld_COMBINED_eSTR_GWASCAT.tab
#ENSGENE=ensgene2.tab
ENSGENE=/storage/resources/dbase/human/gene_annotations/ensgene_to_refseq.tab
MINCAUSAL=0.1

SNPVCF=/storage/resources/datasets/gtex/59533/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v7.p2.c1.GRU/GenotypeFiles/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCaller.vcf.gz
STRVCF=/storage/szfeupe/Runs/650GTEx_estr/Filter_Merged_STRs_All_Samples_New.vcf.gz
CAUSAL=/storage/mgymrek/gtex-estrs-paper/results/eSTR-calling/SuppTable_ALLCAUSAL.csv
