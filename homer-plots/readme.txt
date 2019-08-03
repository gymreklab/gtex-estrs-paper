make_tag_dirs.sh # make tag directories for each BAM file
get_strsets.sh # get BED files with STR lists
get_strsets_GC.sh # get BED files with STR lists for GC analysis
make_composite_nucl.sh # make composite plots for nucleosome positioning
make_composite_nucl_GC.sh # make composite plots for nucleosome positioning for GC analysis
make_composite_RNAPolII_GC.sh # make composite plots for RNAPolII for GC analysis

# Test DNAseI HS GM12878 FM-eSTRs vs. non-eSTRs in 500bp
BAM=/storage/mgymrek/gtex-estrs/revision/homer-plots/encodedata/DNAse_GM12878.bam
cat /storage/mgymrek/gtex-estrs/revision/figures/SuppTable_CAVIAR.tsv | awk -F "\t" '($8>0.3)' | awk '{print $1 "\t" $2-500 "\t" $2+500}' > /storage/mgymrek/gtex-estrs/revision/misc/fmestrs_pm500.bed
cat /storage/mgymrek/gtex-estrs/revision/misc/all_analyzed_strs.tab | cut -f 1,2 | uniq | grep -v chrom | awk '{print $1 "\t" $2-500 "\t" $2+500}' | intersectBed -a stdin -b /storage/mgymrek/gtex-estrs/revision/misc/fmestrs_pm500.bed -v > /storage/mgymrek/gtex-estrs/revision/misc/allestrs_pm500.bed
multiBamCov -bed /storage/mgymrek/gtex-estrs/revision/misc/allestrs_pm500.bed -bams /storage/mgymrek/gtex-estrs/revision/homer-plots/encodedata/DNAse_GM12878.bam > /storage/mgymrek/gtex-estrs/revision/misc/allestrs_pm500_dnase.bed
multiBamCov -bed /storage/mgymrek/gtex-estrs/revision/misc/fmestrs_pm500.bed -bams /storage/mgymrek/gtex-estrs/revision/homer-plots/encodedata/DNAse_GM12878.bam > /storage/mgymrek/gtex-estrs/revision/misc/fmestrs_pm500_dnase.bed
