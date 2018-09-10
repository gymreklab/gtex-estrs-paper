#!/bin/bash

# Get GWAS overlaps
./get_gwas_overlaps.sh

# Filter gwas tables
./filter_gwas.sh GWASCAT

# VLDLR
./run_lz.sh 9 2622147 WholeBlood ENSG00000147852.11 100000
./run_lz_gwas.sh /storage/mgymrek/gtex/gwas/summarystats/hg_Astle_epacs.tab 9_2622147__WholeBlood_ld.txt rs369552432 rs2219143 VLDLR_HGB 2622147
./get_coloc_files_astle_bloodtraits.sh WholeBlood ENSG00000147852.11 9 /storage/mgymrek/gtex/gwas/summarystats/hgb_N172925_wide_form.tsv.gz VLDLR
Rscript --vanilla ./run_coloc_quant.R VLDLR_coloc.txt 172925

# MED19 for SCZ
./run_lz.sh 11 57523883 Adipose-Subcutaneous ENSG00000156603.10 100000
./run_lz_gwas.sh /storage/mgymrek/gtex/gwas/summarystats/pgc_scz2_epacs.tab 11_57523883__Adipose-Subcutaneous_ld.txt rs9420 chr11:57523882 MED19_SCZ2 57523883

# LCAT for HDL
./run_lz.sh 16 68014740 WholeBlood ENSG00000213398.3 100000 
./run_lz_gwas.sh /storage/mgymrek/gtex/gwas/summarystats/pgc_scz2_epacs.tab 16_68014740__WholeBlood_ld.txt rs56303487 chr16:68014739 LCAT_SCZ2 68014740

# DCLK3
./run_lz.sh 3 36835922 Skin-NotSunExposed ENSG00000163673.6 100000
./run_lz_gwas.sh /storage/mgymrek/gtex/gwas/summarystats/pgc_scz2_epacs.tab 3_36835922__Skin-NotSunExposed_ld.txt rs3732386 chr3:36835921 DCLK3_SCZ2 36835922

