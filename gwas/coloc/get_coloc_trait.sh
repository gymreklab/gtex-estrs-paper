#!/bin/bash

source params.sh
PREFIX=$1
TYPE=$2
N=$3

CANDFILE=$OUTDIR/loci/${PREFIX}_candidates.tab
SUMSTATS=$OUTDIR/${PREFIX}_coloc.bed.gz

# For each candidate, get coloc files and perform coloc
while IFS='' read -r line || [[ -n "$line" ]]; do
    gene=$(echo $line | cut -f 1 -d' ')
    tissue=$(echo $line | cut -f 2 -d' ')
    rsid=$(echo $line | cut -f 3 -d' ')
    genename=$(echo $line | cut -f 4 -d' ')

    gtexfile=${OUTDIR}/tmp/${PREFIX}_${gene}_${tissue}_gtex.txt
    gwasfile=${OUTDIR}/tmp/${PREFIX}_${gene}_${tissue}_gwas.txt
    outfile=${OUTDIR}/tmp/${PREFIX}_${gene}_${tissue}_combined.txt
    resfile=${OUTDIR}/${PREFIX}_${genename}_results.txt

    # Get GTEx
    cat /storage/szfeupe/Runs/650GTEx_estr/Analysis_by_Tissue/${tissue}/SNP_Analysis/Lin_Reg_Out | \
	grep -w ${gene} | cut -f 5,10,11,13 | sort -k1,1g > ${gtexfile}
    # Get GWAS
    chrom=$(cat /storage/szfeupe/Runs/650GTEx_estr/Analysis_by_Tissue/${tissue}/SNP_Analysis/Lin_Reg_Out | \
	grep -w ${gene} | cut -f 3 | sed 's/chr//' | head -n 1)
    start=$(cat $gtexfile | datamash min 1)
    end=$(cat $gtexfile | datamash max 1)
    tabix ${SUMSTATS} ${chrom}:${start}-${end} | cut -f 1,3 --complement > ${gwasfile}
    # Combine
    echo "snppos,gwas.beta,gwas.varbeta,gwas.p,gwas.maf,gtex.beta,gtex.varbeta,gtex.p" | sed 's/,/ /g' > ${outfile}
    join --nocheck-order $gwasfile $gtexfile | \
	awk -F" " '{print $1 " " $2 " " $3**2 " " $4 " " $5 " " $6 " " $7**2 " " $8}' >> $outfile
    # Perform coloc
    if [[ x"$TYPE" == x"quant" ]]
    then
	Rscript --vanilla ./run_coloc_quant.R $outfile $N > $resfile
    elif [[ x"$TYPE" == x"cc" ]]
    then
	Rscript --vanilla ./run_coloc_cc.R $outfile $N > $resfile
    else
	echo "invalid type specified $TYPE"
	exit 1
    fi
    # print results to screen
    echo $genename $rsid $tissue $(cat $resfile | head -n 3 | tail -n 1 | awk -F" " '{print $NF}' | sed 's/"//')
done < "$CANDFILE"

