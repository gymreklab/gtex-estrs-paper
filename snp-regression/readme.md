##### Upload necessary data to AWS: #####

1. Expression datasets for each tissue
```
TISSUES="Adipose-Subcutaneous Adipose-Visceral Artery-Aorta Artery-Tibial Brain-Caudate Brain-Cerebellum Cells-Transformedfibroblasts Esophagus-Mucosa Esophagus-Muscularis Heart-LeftVentricle Lung Muscle-Skeletal Nerve-Tibial Skin-NotSunExposed Skin-SunExposed Thyroid WholeBlood"
for tissue in $TISSUES
do
	aws s3 cp /storage/szfeupe/Runs/650GTEx_estr/Analysis_by_Tissue/Review_Rerun/${tissue}/Corr_Expr.csv s3://gtex-estr/Corr_Expr_${tissue}.csv
done
```

2. Expression annotation file
```
aws s3 cp /storage/resources/dbase/human/hg19/gencode_gene_annotations_hg19.csv s3://gtex-estr/gencode_gene_annotations_hg19.csv
```

3. SNP genotype matrices
```
for chrom in $(seq 1 22)
do
	aws s3 cp /storage/mgymrek/gtex-estrs/revision/genotypes/GTExNormalizedSNPGenotypes_chr${chrom}.table.gz s3://gtex-estr/snp_gts_chr${chrom}.tab.gz
	aws s3 cp /storage/mgymrek/gtex-estrs/revision/genotypes/GTExNormalizedSNPGenotypes_chr${chrom}.table.gz.tbi s3://gtex-estr/snp_gts_chr${chrom}.tab.gz.tbi
done
```

##### Make Dockerfile to run regression analysis (1 job/chrom) #####

```
docker build -t gymreklab/gtex-estrs-snpreg .
docker push gymreklab/gtex-estrs-snpreg
```

# Test
```
docker run -it \
       --env AWS_SECRET_ACCESS_KEY=$AWS_SECRET_ACCESS_KEY \
       --env AWS_ACCESS_KEY_ID=$AWS_ACCESS_KEY_ID \
       gymreklab/gtex-estrs-snpreg \
       Nerve-Tibial 1
```

##### Set up AWS Batch environment #####

```
aws batch create-compute-environment \
    --compute-environment-name small \
    --type MANAGED \
    --state ENABLED \
    --compute-resources file://batch-small.json \
    --service-role arn:aws:iam::369425333806:role/service-role/AWSBatchServiceRole

aws batch create-job-queue \
    --job-queue-name gtex-small \
    --state ENABLED \
    --priority 100 \
    --compute-environment-order order=1,computeEnvironment=small

aws batch register-job-definition \
    --job-definition-name gtex-snpreg-job \
    --type container \
    --container-properties file://gtex-snpreg-job.json
```

##### Run AWS jobs #####

```
TISSUES="Adipose-Subcutaneous Adipose-Visceral Artery-Aorta Artery-Tibial Brain-Caudate Brain-Cerebellum Cells-Transformedfibroblasts Esophagus-Mucosa Esophagus-Muscularis Heart-LeftVentricle Lung Muscle-Skeletal Nerve-Tibial Skin-NotSunExposed Skin-SunExposed Thyroid WholeBlood"

for tissue in $TISSUES
do
  #for chrom in $(seq 13 22)
  #for chrom in 3
  do
  cmd="aws batch submit-job \
      --job-name ${tissue}-${chrom} \
      --job-queue gtex-small \
      --job-definition gtex-snpreg-job:3 \
      --container-overrides 'command=[\"${tissue}\",\"${chrom}\"]'"
  sh -c "${cmd}"
  done 
done
```