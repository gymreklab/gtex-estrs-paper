# Build docker
```
docker build -t gymreklab/gtex-estrs-strreg .
docker push gymreklab/gtex-estrs-strreg
```

# Set up job
```
aws batch register-job-definition \
    --job-definition-name gtex-strreg-job \
    --type container \
    --container-properties file://gtex-strreg-job.json
```

# Test AWS
```
tissue=Brain-Cerebellum
chrom=3
  cmd="aws batch submit-job \
      --job-name ${tissue}-${chrom} \
      --job-queue gtex-small \
      --job-definition gtex-strreg-job:1 \
      --container-overrides 'command=[\"${tissue}\",\"${chrom}\"]'"
  sh -c "${cmd}"
```

# Run all AWS
```
TISSUES="Adipose-Subcutaneous Adipose-Visceral Artery-Aorta Artery-Tibial Brain-Caudate Brain-Cerebellum Cells-Transformedfibroblasts Esophagus-Mucosa Esophagus-Muscularis Heart-LeftVentricle Lung Muscle-Skeletal Nerve-Tibial Skin-NotSunExposed Skin-SunExposed Thyroid WholeBlood"

for tissue in $TISSUES
do
	for chrom in $(seq 1 22)
  do
  cmd="aws batch submit-job \
      --job-name ${tissue}-${chrom} \
      --job-queue gtex-small \
      --job-definition gtex-strreg-job:1 \
      --container-overrides 'command=[\"${tissue}\",\"${chrom}\"]'"
  sh -c "${cmd}"
  done 
done
```