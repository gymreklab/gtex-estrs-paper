#!/bin/bash

# Get summary stats ready for coloc
./make_coloc_summstats.sh

# Get list of gene/tissues for each trait and perform coloc
./get_candidates.sh Height height
./get_coloc_trait.sh height quant 695647 > ${OUTDIR}/height_coloc_combined.tab
./summarize_coloc.sh height > ${OUTDIR}/height_coloc_summary.tab

./get_candidates.sh Schizophrenia scz_pgc
./get_coloc_trait.sh scz_pgc cc 0.33 > ${OUTDIR}/scz_pgc_coloc_combined.tab
./summarize_coloc.sh scz_pgc > ${OUTDIR}/scz_pgc_coloc_summary.tab

./get_candidates.sh "Hemoglobin concentration" hgb_N172925
./get_coloc_trait.sh hgb_N172925 quant 172925 > ${OUTDIR}/hgb_coloc_combined.tab
./summarize_coloc.sh hgb_N172925 > ${OUTDIR}/hgb_coloc_summary.tab

./get_candidates.sh "Basophil percentage of granulocytes" baso_p_gran_N170223
./get_coloc_trait.sh baso_p_gran_N170223 quant 171996 > ${OUTDIR}/baso_p_gran_combined.tab

./get_candidates.sh "Red blood cell count" rbc_N172952
./get_coloc_trait.sh rbc_N172952 quant 172952 > ${OUTDIR}/rbc_combined.tab

#./get_candidates.sh "" hdl

# Get LD for the ones we want to plot
./get_ld.sh hgb_N172925 ENSG00000147852.11 WholeBlood 9:2622147
./get_ld.sh scz_pgc ENSG00000156603.10 Adipose-Subcutaneous 11:57523883
./get_ld.sh scz_pgc ENSG00000126464.9 Thyroid 19:50110619
./get_ld.sh scz_pgc ENSG00000143641.8 Artery-Aorta 1:230264005
./get_ld.sh height ENSG00000163933.5 Artery-Aorta 3:53128363
./get_ld.sh height ENSG00000115556.9 Muscle-Skeletal 2:219482301
./get_ld.sh scz_pgc ENSG00000112893.5 WholeBlood 5:109038017
./get_ld.sh scz_pgc ENSG00000213996.8 Adipose-Visceral 19:19424949
