#!/bin/bash

TISSUES="Adipose-Subcutaneous Adipose-Visceral Artery-Aorta Artery-Tibial Brain-Caudate Brain-Cerebellum Cells-Transformedfibroblasts Esophagus-Mucosa Esophagus-Muscularis Heart-LeftVentricle Lung Muscle-Skeletal Nerve-Tibial Skin-NotSunExposed Skin-SunExposed Thyroid WholeBlood"

# Number genes tested per and number of STR*gene tests per tissue
for tissue in $TISSUES
do
    tfile=/storage/mgymrek/gtex-estrs/revision/strreg/${tissue}_strreg.tab
    echo $tissue $(cat $tfile | grep -v gene | wc -l) $(cat $tfile | grep -v gene | cut -f 1 | sort | uniq | wc -l)
done

exit 1

# Number of STRs tested per gene per tissue
for tissue in $TISSUES
do
    cat /storage/mgymrek/gtex-estrs/revision/strreg/${tissue}_strreg.tab | \
	grep -v gene | \
	cut -f 1 | uniq -c | awk '{print $1 "\t" $2}'
done | datamash mean 1

