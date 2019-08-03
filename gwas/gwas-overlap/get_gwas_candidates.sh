#!/bin/bash

SNPVCF=/storage/resources/datasets/gtex/59533/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v7.p2.c1.GRU/GenotypeFiles/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCaller.vcf.gz
STRVCF=/storage/szfeupe/Runs/650GTEx_estr/Filter_Merged_STRs_All_Samples_New.vcf.gz \

for trait in GlaucomaUKBB #SCZPGC HeightYengo IBDHuang IntelligenceSavageJensen AlzheimersIGAP
do
    filename=/storage/mgymrek/gtex-estrs/revision/coloc/gwashits/${trait}.tab
    ./GetGWASOverlap.py \
	--estrs /storage/mgymrek/gtex-estrs/revision/figures/SuppTable_CAVIAR.tsv \
	--gwas ${filename} --prefix ${trait} \
	--window 1000000 \
	--caviar 0.3 \
	--outdir /storage/mgymrek/gtex-estrs/revision/coloc/candidates/ \
	--str ${STRVCF} \
	--snp ${SNPVCF}
done
