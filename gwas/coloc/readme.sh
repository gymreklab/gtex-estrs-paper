# See make_coloc_summstats.sh to get summ stats in necessary format

# Run full coloc analysis - with hits from each study
./run_trait.sh SCZPGC cc 0.33
./run_trait.sh HeightYengo quant 695647
./run_trait.sh IBDHuang cc 0.33
./run_trait.sh IntelligenceSavageJensen quant 269720
./run_trait.sh AlzheimersIGAP cc 0.33

# Run full coloc analysis - with NHGRI hits
./run_trait_nhgri.sh SCZPGC cc 0.33 
./run_trait_nhgri.sh HeightYengo quant 695647
./run_trait_nhgri.sh IBDHuang cc 0.33
./run_trait_nhgri.sh IntelligenceSavageJensen quant 269720
./run_trait_nhgri.sh AlzheimersIGAP cc 0.33

# Get coloc files for plotting GWAS hits. See make_files.txt for details
./get_coloc_trait.sh SIGLEC14 NA NA

OUTDIR=/storage/mgymrek/gtex-estrs/revision/coloc/nhgri/results/
for trait in IBDHuang IntelligenceSavageJensen SCZPGC HeightYengo
do
    ./summarize_coloc_nhgri.sh ${trait} > ${OUTDIR}/${trait}_coloc_results.tab
done

for trait in IBDHuang IntelligenceSavageJensen SCZPGC HeightYengo
do
    cat ${OUTDIR}/${trait}_coloc_results.tab | grep -v gene | \
	awk -F"\t" -v"trait=$trait" '{print $1":"$5 "\t" $2 "\t" $3 "\t" $4 "\t" $6 "\t" $7"\t" $8 "\t" $13 "\t" trait}'
done > /storage/mgymrek/gtex-estrs/revision/figures/SuppTable_Coloc.tsv
