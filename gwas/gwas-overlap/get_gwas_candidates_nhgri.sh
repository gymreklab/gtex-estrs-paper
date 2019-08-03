#!/bin/bash

SNPVCF=/storage/resources/datasets/gtex/59533/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v7.p2.c1.GRU/GenotypeFiles/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCaller.vcf.gz
STRVCF=/storage/szfeupe/Runs/650GTEx_estr/Filter_Merged_STRs_All_Samples_New.vcf.gz \
GWAS=/storage/mgymrek/gtex/gwas/summarystats/ucsc_gwas_catalog_072419_v2.tab

# Schizophrenia
grep "^chrom\|Schizophrenia" ${GWAS} > /storage/mgymrek/gtex/gwas/summarystats/GWASCatalog_Schizophrenia.tab
./GetGWASOverlap.py \
    --estrs /storage/mgymrek/gtex-estrs/revision/figures/SuppTable_CAVIAR.tsv \
    --gwas /storage/mgymrek/gtex/gwas/summarystats/GWASCatalog_Schizophrenia.tab \
    --prefix SCZPGC_nhgri \
    --window 1000000 \
    --caviar 0.3 \
    --outdir /storage/mgymrek/gtex-estrs/revision/coloc/candidates/ \
    --str ${STRVCF} \
    --snp ${SNPVCF}

# Height
grep "^chrom\|Height" ${GWAS} > /storage/mgymrek/gtex/gwas/summarystats/GWASCatalog_Height.tab
./GetGWASOverlap.py \
    --estrs /storage/mgymrek/gtex-estrs/revision/figures/SuppTable_CAVIAR.tsv \
    --gwas /storage/mgymrek/gtex/gwas/summarystats/GWASCatalog_Height.tab \
    --prefix HeightYengo_nhgri \
    --window 1000000 \
    --caviar 0.3 \
    --outdir /storage/mgymrek/gtex-estrs/revision/coloc/candidates/ \
    --str ${STRVCF} \
    --snp ${SNPVCF}

# IBD
grep "^chrom\|Ulcerative\|Inflammatory\|Crohn" ${GWAS} | grep -v psoriasis > /storage/mgymrek/gtex/gwas/summarystats/GWASCatalog_IBD.tab
./GetGWASOverlap.py \
    --estrs /storage/mgymrek/gtex-estrs/revision/figures/SuppTable_CAVIAR.tsv \
    --gwas /storage/mgymrek/gtex/gwas/summarystats/GWASCatalog_IBD.tab \
    --prefix IBDHuang_nhgri \
    --window 1000000 \
    --caviar 0.3 \
    --outdir /storage/mgymrek/gtex-estrs/revision/coloc/candidates/ \
    --str ${STRVCF} \
    --snp ${SNPVCF}

# Intelligence
grep "^chrom\|Intelligence" ${GWAS} > /storage/mgymrek/gtex/gwas/summarystats/GWASCatalog_Intelligence.tab
./GetGWASOverlap.py \
    --estrs /storage/mgymrek/gtex-estrs/revision/figures/SuppTable_CAVIAR.tsv \
    --gwas /storage/mgymrek/gtex/gwas/summarystats/GWASCatalog_Intelligence.tab \
    --prefix IntelligenceSavageJensen_nhgri \
    --window 1000000 \
    --caviar 0.3 \
    --outdir /storage/mgymrek/gtex-estrs/revision/coloc/candidates/ \
    --str ${STRVCF} \
    --snp ${SNPVCF}

# Alzheimer
grep "^chrom\|Alzheimer" ${GWAS} > /storage/mgymrek/gtex/gwas/summarystats/GWASCatalog_Alzheimer.tab
./GetGWASOverlap.py \
    --estrs /storage/mgymrek/gtex-estrs/revision/figures/SuppTable_CAVIAR.tsv \
    --gwas /storage/mgymrek/gtex/gwas/summarystats/GWASCatalog_Intelligence.tab \
    --prefix AlzheimersIGAP_nhgri \
    --window 1000000 \
    --caviar 0.3 \
    --outdir /storage/mgymrek/gtex-estrs/revision/coloc/candidates/ \
    --str ${STRVCF} \
    --snp ${SNPVCF}
