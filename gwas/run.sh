#!/bin/bash

source params.sh

GWASBED=$1 # chrom, start, rsid, locusname, pval
PREFIX=$2

##### Get all STRs within window kb of a GWAS catalog SNP
echo "Get all nearby STRs that are also eSTRs...."
cat ${CAUSAL} | sed 's/,/\t/g' | grep -v gene | awk '{print $2 "\t" $3-1 "\t" $3+1}' | sed 's/chr//' > ${OUTDIR}/tmp/all_estrs.bed
cat ${GWASBED} | grep -v chrom | awk '{print $1 "\t" $2 "\t" $2+1 "\t" $0}' > ${OUTDIR}/tmp/${PREFIX}.bed
cat ${HIPREF} | sed 's/chr//' | \
    intersectBed -a stdin -b ${OUTDIR}/tmp/all_estrs.bed -wa | \
    awk -v"window=$WINDOW" '{print $1 "\t" ($2<window?0:$2-window) "\t" $3+window "\t" $0}' | \
    intersectBed -a stdin -b ${OUTDIR}/tmp/${PREFIX}.bed -wa -wb | \
    cut -f 1-3 --complement  > ${OUTDIR}/tmp/str_gwas_overlap_${PREFIX}.bed
cat ${OUTDIR}/tmp/str_gwas_overlap_${PREFIX}.bed | cut -f 1,2,7 | sort | uniq > ${OUTDIR}/tmp/str_gwas_coords_${PREFIX}.tab

##### Get SNP-STR LD
echo "Get SNP-STR LD..."
~/workspace/ssc-imputation/snpstr-ld/snp_str_ld_calculator.py \
    --return-filtered-snps \
    --str-vcf ${STRVCF} \
    --snp-vcf ${SNPVCF} \
    --loci-file ${OUTDIR}/tmp/str_gwas_coords_${PREFIX}.tab \
    --use-info-start --mincount 3 --usefilter --use-gb > ${OUTDIR}/tmp/str_gwas_ld_hg19_${PREFIX}.tab

##### Combine
echo "Combine..."
cat ${OUTDIR}/tmp/str_gwas_ld_hg19_${PREFIX}.tab | grep -v locus1 | sed 's/:/\t/g' | \
    awk '{print $1 "\t" $2 "\t" $2+1 "\t" $0}' | intersectBed -a stdin -b ${OUTDIR}/tmp/str_gwas_overlap_${PREFIX}.bed -wa -wb | \
    awk '(($2==$15) && ($7==$20))' | cut -f 4,5,12,17,18,20,24- > ${OUTDIR}/str_gwas_ld_COMBINED_${PREFIX}.tab

# ,chrom,str.start,str.end,gene,gene.name,dist.tss,motif,score,beta,pval,qval,num.e,causal,tissue_info
##### Combine with eSTRs
cat ${OUTDIR}/str_gwas_ld_COMBINED_${PREFIX}.tab | awk '{print $1 "\t" $2 "\t" $2+1 "\t" $0}' | sed 's/chr//' > ${OUTDIR}/tmp/del1${PREFIX}
cat ${CAUSAL} | sed 's/,/\t/g' | \
    grep -v chrom | \
    awk '{print $2 "\t" $3-1 "\t" $4-1 "\t" $6 "\t" $9 "\t" $10 "\t" $15}' | sed 's/chr//' > ${OUTDIR}/tmp/del2${PREFIX}
echo "chrom,start,LD,period,motif,snp.start,rsid,locus,gene.name,score,beta,tissue_info" | sed 's/,/\t/g' \
    > ${OUTDIR}/str_gwas_ld_COMBINED_eSTR_${PREFIX}.tab 
intersectBed -a ${OUTDIR}/tmp/del1${PREFIX} -b ${OUTDIR}/tmp/del2${PREFIX} -wa -wb  | cut -f 4-11,15- | uniq \
    >> ${OUTDIR}/str_gwas_ld_COMBINED_eSTR_${PREFIX}.tab 


