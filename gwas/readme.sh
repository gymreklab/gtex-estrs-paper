#!/bin/bash

# Get GWAS overlaps
./get_gwas_overlaps.sh

# VLDLR
./run_lz.sh 9 2622147 WholeBlood ENSG00000147852.11 100000
./run_lz_gwas.sh /storage/mgymrek/gtex/gwas/summarystats/hg_Astle_epacs.tab 9_2622147__WholeBlood_ld.txt rs369552432 rs2219143 VLDLR_HGB 2622147

# MED19 for SCZ
./run_lz.sh 11 57523883 Adipose-Subcutaneous ENSG00000156603.10 100000
./run_lz_gwas.sh /storage/mgymrek/gtex/gwas/summarystats/pgc_scz2_epacs.tab 11_57523883__Adipose-Subcutaneous_ld.txt rs9420 chr11:57523882 MED19_SCZ2 57523883
