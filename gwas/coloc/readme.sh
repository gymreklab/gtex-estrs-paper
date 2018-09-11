#!/bin/bash

# Get summary stats ready for coloc
./make_coloc_summstats.sh

# Get list of gene/tissues for each trait and perform coloc
./get_candidates.sh Schizophrenia scz_pgc
./get_coloc_trait.sh scz_pgc cc 0.33 > ${OUTDIR}/scz_pgc_coloc_combined.tab

./get_candidates.sh "Hemoglobin concentration" hgb_N172925
./get_coloc_trait.sh hgb_N172925 quant 172925 > ${OUTDIR}/hgb_coloc_combined.tab

./get_candidates.sh "Basophil percentage of granulocytes" baso_p_gran_N170223
./get_coloc_trait.sh baso_p_gran_N170223 quant 171996 > ${OUTDIR}/baso_p_gran_combined.tab

./get_candidates.sh "Red blood cell count" rbc_N172952
./get_coloc_trait.sh rbc_N172952 quant 172952 > ${OUTDIR}/rbc_combined.tab

# Get LD for the ones we want to plot
./get_ld.sh hgb_N172925 ENSG00000147852.11 WholeBlood 9:2622147
