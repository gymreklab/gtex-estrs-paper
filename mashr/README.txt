### How to run ###

# Prepare input data
./gatherData.sh

# Run - STRs
nohup R CMD BATCH '--args runval="strs"' runMashr.R /storage/mgymrek/gtex-estrs/revision/mashr/output-strs/mashr.strs.log &

# Run - SNPs
# Learn model on chr1
nohup R CMD BATCH '--args runval="snps-bychrom/chr1"' runMashr.R /storage/mgymrek/gtex-estrs/revision/mashr/output-snps-bychrom/chr1/mashr.model.log &
# Run all chroms separately, using chr1 model, in chunks of 10000
nohup ./run_mashr_snps.sh &

# Collate by chrom
./collate-chroms.sh /storage/mgymrek/gtex-estrs/revision/mashr/ strs
./collate-chroms.sh /storage/mgymrek/gtex-estrs/revision/mashr/ snps

# Compute Z-scores
./compute-mashR-Z.py /storage/mgymrek/gtex-estrs/revision/mashr/ strs
./compute-mashR-Z.py /storage/mgymrek/gtex-estrs/revision/mashr/ snps

# Compute significant eSTRs/eSNPs
./compute-mashR-sig.py /storage/mgymrek/gtex-estrs/revision/mashr/ strs
./compute-mashR-sig.py /storage/mgymrek/gtex-estrs/revision/mashr/ snps

# Identify tissue-specific eSTRs
./find-tissue-specific.py /storage/mgymrek/gtex-estrs/revision/mashr/output-strs/zscores.tsv /storage/mgymrek/gtex-estrs/revision/mashr/output-strs/tissue-specific
