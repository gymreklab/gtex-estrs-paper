#!/bin/bash

# Hemoglobin Astle et al
zcat /storage/mgymrek/gtex/gwas/summarystats/hgb_N172925_wide_form.tsv.gz | awk '{print $4 "\t" $5 "\t" $5 "\t" $3 "\tNS\tAC\tCALLRATE\tMAF\t" $13 "\tSCORE\tNCASE\tNCTRL\tAFCASE\tAFCTRL"}' | grep -v CHR >> /storage/mgymrek/gtex/gwas/summarystats/hgb_Astle_epacs.tab

# SCZ
#head -n 1 /storage/mgymrek/gtex/gwas/summarystats/hgb_Astle_epacs.tab > pgc_scz2_epacs.tab∑
#zcat SCZvsCONT.sumstats.gz | awk '{print $1 "\t" $3 "\t" $3 "\t" $2 "\tNS\tAC\tCALLRATE\tMAF\t" $11 "\tSCORE\tNCASE\tNCTRL\tAFCASE\tAFCTRL" }' | grep -v CHR >> pgc_scz2_epacs.tab

head -n 1 /storage/mgymrek/gtex/gwas/summarystats/hgb_Astle_epacs.tab > /storage/mgymrek/gtex/gwas/summarystats/pgc_scz2_epacs.tab
zcat /storage/mgymrek/gtex/gwas/summarystats/ckqny.scz2snpres.gz | grep -v hg19 | \
    awk '{print $1 "\t" $5 "\t" $5 "\t" $2 "\tNS\tAC\tCALLRATE\tMAF\t" $9 "\tSCORE\tNCASE\tNCTRL\tAFCASE\tAFCTRL" }' | sed 's/chr//' >> /storage/mgymrek/gtex/gwas/summarystats/pgc_scz2_epacs.tab
