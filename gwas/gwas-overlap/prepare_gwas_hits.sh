#!/bin/bash

# Get GWAS hits in format chrom start rsid name
# Also put IntelligenceSavageJensen.tab AlzheimersIGAP.tab manually

OUTDIR=/storage/mgymrek/gtex-estrs/revision/coloc/gwashits

# SCZ
echo "chrom,start,rsid,name" | sed 's/,/\t/g' > ${OUTDIR}/SCZPGC.tab
cat /storage/mgymrek/gtex/gwas/summarystats/PGC_SCZ2_hits_v2.txt | grep -v name | \
    awk '{print $1 "\t" $2 "\t" $3 "\t" "SCZ-PGC"}' >> ${OUTDIR}/SCZPGC.tab

# Height
echo "chrom,start,rsid,name" | sed 's/,/\t/g' > ${OUTDIR}/HeightYengo.tab
cat /storage/mgymrek/gtex/gwas/summarystats/Height_Yengo_2018.txt | grep -v name | \
    awk '{print $1 "\t" $2 "\t" $3 "\t" "Height-Yengo"}' >> ${OUTDIR}/HeightYengo.tab

# IBD
cat /storage/mgymrek/gtex/gwas/summarystats/Huang_etal_IBD_GWAShits.tab.txt  | \
    cut -f 1-4 > ${OUTDIR}/IBDHuang.tab
