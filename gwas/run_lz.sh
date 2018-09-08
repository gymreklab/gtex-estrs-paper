#!/bin/bash

# Example: ./run_lz.sh 10 134000824 Cells-Transformedfibroblasts ENSG00000151640.8 100000 
# Example: ./run_lz.sh 9 2622147 WholeBlood ENSG00000147852.11 100000

source params.sh

CHROM=$1
START=$2
TISSUE=$3
GENE=$4
WINDOW=$5

PFILE=${CHROM}_${START}_${GEHE}_${TISSUE}_pvals.txt
SFILE=${CHROM}_${START}_${GEHE}_${TISSUE}_snpstr.txt
LFILE=${CHROM}_${START}_${GEHE}_${TISSUE}_ld.txt

# Pull out data in epacts format
echo "Pulling out pvalues..."
# The chrom, start, end, marker ID, and p-value columns must all be present. The file must be tab-delimited.
echo "#CHROM,BEGIN,END,MARKER_ID,NS,AC,CALLRATE,MAF,PVALUE,SCORE,N.CASE,N.CTRL,AF.CASE,AF.CTRL" | sed 's/,/\t/g'> ${PFILE}

cat /storage/szfeupe/Runs/650GTEx_estr/Analysis_by_Tissue/${TISSUE}/Lin_Reg_Out | \
    grep -w ${START} | grep -w ${GENE} | sed 's/chr//' | \
    awk '{print $2 "\t" $4 "\t" $4 "\t" $2":"$4 "\t" "NS\tAC\tCALLRATE\tMAF\t" $12 "\t" "SCORE\tNCASE\tNCTRL\tAFCASE\tAFCTRL"}' >> ${PFILE}
cat /storage/szfeupe/Runs/650GTEx_estr/Analysis_by_Tissue/${TISSUE}/SNP_Analysis/Lin_Reg_Out | cut -f 1 --complement | \
    grep -w ${GENE} | sed 's/chr//' | \
    awk '{print $2 "\t" $4 "\t" $4 "\t" $2":"$4  "\t" "NS\tAC\tCALLRATE\tMAF\t" $12 "\t" "SCORE\tNCASE\tNCTRL\tAFCASE\tAFCTRL"}' |
    sed 's/STR_/SNP_/' >> ${PFILE}
    
# Compute LD
echo "Computing LD..."
cat ${PFILE} | grep -v "STR" | grep -v "CHROM" | cut -f 2 | \
    awk -v"chrom=$CHROM" -v"start=$START" '{print chrom "\t" start "\t" $0}' > ${SFILE}
echo "snp1,snp2,dprime,rsquare" | sed 's/,/\t/g' > ${LFILE}
echo "chr${CHROM}:$((START-1)),chr${CHROM}:$((START-1)),NA,1" | sed 's/,/\t/g' >> ${LFILE}
#echo "STR_${START},STR_${START},NA,1" | sed 's/,/\t/g' >> ${LFILE}
/home/mgymrek/workspace/ssc-imputation/snpstr-ld/snp_str_ld_calculator.py \
    --str-vcf ${STRVCF} --snp-vcf ${SNPVCF} --loci-file ${SFILE} \
    --use-info-start --mincount 3 --usefilter --use-gb | \
    grep -v locus1 | sed 's/:/\t/g' | \
    awk -v"chrom=$CHROM" '{print "chr"chrom":"$4 "\t" "chr"chrom":"$2-1 "\tNA\t" $9}' >> ${LFILE}
    
# Run locus zoom on expression data
echo "Run LZ..."
/storage/resources/source/locuszoom/bin/locuszoom \
    --epacts ${PFILE} --refsnp chr${CHROM}:$((START-1)) \
    --chr ${CHROM} --start $((START-$WINDOW)) --end $((START+$WINDOW)) \
    --ld ${LFILE} --build hg19 \
    --prefix lz/${GENE}_${CHROM}_${START}


