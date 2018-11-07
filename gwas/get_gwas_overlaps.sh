#!/bin/bash

# GWAS catalog
./run.sh /storage/mgymrek/gtex/gwas/summarystats/ucsc_gwas_catalog_090818_v2.tab GWASCAT3

# SCZ
./run.sh /storage/mgymrek/gtex/gwas/summarystats/PGC_SCZ2_hits_v2.txt SCZ2

# Height
./run.sh /storage/mgymrek/gtex/gwas/summarystats/Height_Yengo_2018.txt Height

# IBD
./run.sh /storage/mgymrek/gtex/gwas/summarystats/Huang_etal_IBD_GWAShits.tab.txt IBD
