#!/bin/bash

# VLDLR
./run_lz.sh 9 2622147 WholeBlood ENSG00000147852.11 100000
./run_lz_gwas.sh /storage/mgymrek/gtex/gwas/summarystats/hg_Astle_epacs.tab 9_2622147__WholeBlood_ld.txt rs369552432 rs2219143 VLDLR_HGB
