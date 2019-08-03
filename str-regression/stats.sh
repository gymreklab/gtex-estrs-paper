#!/bin/bash

# Number of STRs tested per gene per tissue
TISSUES="Adipose-Subcutaneous Adipose-Visceral Artery-Aorta Artery-Tibial Brain-Caudate Brain-Cerebellum Cells-Transformedfibroblasts Esophagus-Mucosa Esophagus-Muscularis Heart-LeftVentricle Lung Muscle-Skeletal Nerve-Tibial Skin-NotSunExposed Skin-SunExposed Thyroid WholeBlood"
for tissue in $TISSUES
do
    cat /storage/mgymrek/gtex-estrs/revision/strreg/${tissue}_strreg.tab | \
	grep -v gene | \
	cut -f 1 | uniq -c | awk '{print $1 "\t" $2}'
done | datamash mean 1
