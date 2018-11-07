#!/bin/bash

# Get summary stats ready for coloc
./make_coloc_summstats.sh

# Get list of gene/tissues for each trait and perform coloc
./get_candidates.sh Height height
./get_coloc_trait.sh height quant 695647 > ${OUTDIR}/height_coloc_combined.tab
./summarize_coloc.sh height Height > ${OUTDIR}/height_coloc_summary.tab
./plot_candidates.sh height
cat /storage/mgymrek/gtex/gwas/summarystats/coloc/height_coloc_summary.tab  | grep -v gene | awk '($12>0.5)' | sort -k6,6g | awk '{print $1":"$4 "\t" $2 "\t" $3 "\t" $5 "\t" $6 "\t" $7 "\t" $12}'

./get_candidates.sh SCZ2 scz_pgc
./get_coloc_trait.sh scz_pgc cc 0.33 > ${OUTDIR}/scz_pgc_coloc_combined.tab
./summarize_coloc.sh scz_pgc SCZ2 > ${OUTDIR}/scz_pgc_coloc_summary.tab
./plot_candidates.sh scz_pgc
cat /storage/mgymrek/gtex/gwas/summarystats/coloc/scz_pgc_coloc_summary.tab  | grep -v gene | awk '($12>0.5)' | sort -k6,6g | awk '{print $1":"$4 "\t" $2 "\t" $3 "\t" $5 "\t" $6"\t" $7 "\t" $12}'

./get_candidates.sh IBD ibd
./get_coloc_trait.sh ibd cc 0.33 > ${OUTDIR}/ibd_coloc_combined.tab
./summarize_coloc.sh ibd IBD > ${OUTDIR}/ibd_coloc_summary.tab
./plot_candidates.sh ibd
cat /storage/mgymrek/gtex/gwas/summarystats/coloc/ibd_coloc_summary.tab  | grep -v gene | awk '($12>0.5)' | sort -k6,6g | awk '{print $1":"$4 "\t" $2 "\t" $3 "\t" $5 "\t" $6 "\t" $7 "\t" $12}'
