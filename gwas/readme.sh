#!/bin/bash

# Get GWAS overlaps
./get_gwas_overlaps.sh

# Filter gwas tables
./filter_gwas.sh GWASCAT

########### Hemoglobin #########
# VLDLR
N_HGB=172925 # See Astle table S2
HGBSTATS=/storage/mgymrek/gtex/gwas/summarystats/hgb_N172925_wide_form.tsv.gz

./run_lz.sh 9 2622147 WholeBlood ENSG00000147852.11 100000
./run_lz_gwas.sh /storage/mgymrek/gtex/gwas/summarystats/hg_Astle_epacs.tab 9_2622147__WholeBlood_ld.txt rs369552432 rs2219143 VLDLR_HGB 2622147


./get_coloc_files_astle_bloodtraits.sh WholeBlood ENSG00000147852.11 9 ${HGBSTATS} VLDLR
Rscript --vanilla ./run_coloc_quant.R VLDLR_coloc.txt $N_HGB

# ESR2
./get_coloc_files_astle_bloodtraits.sh Skin-SunExposed ENSG00000140009.14 14 ${HGBSTATS} ESR2
Rscript --vanilla ./run_coloc_quant.R ESR2_coloc.txt $N_HGB

# SLC36A1
BSPSTATS=/storage/mgymrek/gtex/gwas/summarystats/baso_p_N171996_wide_form.tsv.gz
N_BSP=171996
./get_coloc_files_astle_bloodtraits.sh Thyroid ENSG00000123643.8 5 ${BSPSTATS} SLC36A1
Rscript --vanilla ./run_coloc_quant.R SLC36A1_coloc.txt $N_BSP

############ Platelet count ###########
# https://www.nature.com/articles/s41588-018-0047-6
PLATELETSTATS=/storage/mgymrek/gtex/gwas/summarystats/BBJ.Plt.autosome.txt.gz
N_PLATELET=108208
./get_coloc_files_japanese.sh Cells-Transformedfibroblasts ENSG00000160703.11 11 ${PLATELETSTATS} NLRX1
Rscript --vanilla ./run_coloc_quant.R NLRX1_coloc.txt $N_PLATELET

##############################

# APH1A for outercanthal width
./run_lz.sh 1 150314982 Cells-Transformedfibroblasts ENSG00000117362.8 100000

# MED19 for SCZ
./run_lz.sh 11 57523883 Adipose-Subcutaneous ENSG00000156603.10 100000
./run_lz_gwas.sh /storage/mgymrek/gtex/gwas/summarystats/pgc_scz2_epacs.tab 11_57523883__Adipose-Subcutaneous_ld.txt rs9420 chr11:57523882 MED19_SCZ2 57523883

# LCAT for HDL
./run_lz.sh 16 68014740 WholeBlood ENSG00000213398.3 100000 
./run_lz_gwas.sh /storage/mgymrek/gtex/gwas/summarystats/pgc_scz2_epacs.tab 16_68014740__WholeBlood_ld.txt rs56303487 chr16:68014739 LCAT_SCZ2 68014740

# DCLK3
./run_lz.sh 3 36835922 Skin-NotSunExposed ENSG00000163673.6 100000
./run_lz_gwas.sh /storage/mgymrek/gtex/gwas/summarystats/pgc_scz2_epacs.tab 3_36835922__Skin-NotSunExposed_ld.txt rs3732386 chr3:36835921 DCLK3_SCZ2 36835922

