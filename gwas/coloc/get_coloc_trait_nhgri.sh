#!/bin/bash

trait=$1
type=$2
n=$3

OUTDIR=/storage/mgymrek/gtex-estrs/revision/coloc/nhgri
CANDFILE=/storage/mgymrek/gtex-estrs/revision/coloc/candidates/${trait}_nhgri_candidates.tab
SUMSTATS=/storage/mgymrek/gtex-estrs/revision/coloc/summstats/${trait}_coloc.bed.gz

# For each candidate, get coloc files and perform coloc
while IFS='' read -r line || [[ -n "$line" ]]; do
    gene=$(echo $line | cut -f 1 -d' ')
    tissue=$(echo $line | cut -f 2 -d' ')
    genename=$(echo $line | cut -f 3 -d' ')

    if [ x"$tissue" == x"tissue" ]; then continue; fi

    gtexfile=${OUTDIR}/tmp/${trait}_${gene}_${tissue}_gtex.txt
    gwasfile=${OUTDIR}/tmp/${trait}_${gene}_${tissue}_gwas.txt
    outfile=${OUTDIR}/tmp/${trait}_${gene}_${tissue}_combined.txt
    resfile=${OUTDIR}/res/${trait}_${genename}_results.txt

    # Get GTEx. Need str.start,beta,beta.se,p.wald
    cat /storage/mgymrek/gtex-estrs/revision/snpreg/${tissue}_snpreg.tab | \
	grep -w ${gene} | cut -f 4,6,7,9 | sort -k1,1g > ${gtexfile}
    # Get coordinates
    chrom=$(cat /storage/mgymrek/gtex-estrs/revision/snpreg/${tissue}_snpreg.tab | \
	grep -w ${gene} | cut -f 2 | sed 's/chr//' | head -n 1)
    start=$(cat $gtexfile | datamash min 1)
    end=$(cat $gtexfile | datamash max 1)
    # Get GWAS
    tabix ${SUMSTATS} ${chrom}:${start}-${end} | cut -f 1,3 --complement > ${gwasfile}
    # Combine
    echo "snppos,gwas.beta,gwas.varbeta,gwas.p,gwas.maf,gtex.beta,gtex.varbeta,gtex.p" | sed 's/,/ /g' > ${outfile}
    join --nocheck-order $gwasfile $gtexfile | \
	awk -F" " '{print $1 " " $2 " " $3**2 " " $4 " " $5 " " $6 " " $7**2 " " $8}' >> $outfile
    # Perform coloc
    if [[ x"$type" == x"quant" ]]
    then
	Rscript --vanilla ./run_coloc_quant.R $outfile $n > $resfile
    elif [[ x"$type" == x"cc" ]]
    then
	Rscript --vanilla ./run_coloc_cc.R $outfile $n > $resfile
    else
	echo "invalid type specified $type"
	exit 1
    fi
    # print results to screen
    echo $genename $tissue $(cat $resfile | head -n 3 | tail -n 1 | awk -F" " '{print $NF}' | sed 's/"//')
done < "$CANDFILE"

