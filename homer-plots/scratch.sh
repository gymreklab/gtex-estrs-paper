# Get per-tissue eSTRs by period
TISSUES="Adipose-Subcutaneous Adipose-Visceral Artery-Aorta Artery-Tibial Brain-Caudate Brain-Cerebellum Cells-Transformedfibroblasts Esophagus-Mucosa Esophagus-Muscularis Heart-LeftVentricle Lung Muscle-Skeletal Nerve-Tibial Skin-NotSunExposed Skin-SunExposed Thyroid WholeBlood"
for tissue in $TISSUES
do
    for period in $(seq 1 6)
    do
	cat ${MASTER}/${tissue}_master.tab | \
	    csvcut -t -c chrom,str.start,str.end,str.motif.forward,caviar.str.score | \
	    csvformat -T | grep -v chrom | \
	    awk -v"period=$period" -v"score=$SCORE" \
	    '(length($4)==period && $5>=score) {print $1 "\t" $2 "\t" $3}' | \
	    sort | uniq > ${OUTDIR}/${tissue}_period${period}.bed
    done
    cat ${MASTER}/${tissue}_master.tab | \
	csvcut -t -c chrom,str.start,str.end,str.motif.forward,caviar.str.score | \
	csvformat -T | grep -v chrom | \
	awk  -v"score=$SCORE" \
	'($5>=score) {print $1 "\t" $2 "\t" $3}' | \
	sort | uniq > ${OUTDIR}/${tissue}_ALL.bed
done 
