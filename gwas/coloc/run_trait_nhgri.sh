#!/bin/bash

OUTDIR=/storage/mgymrek/gtex-estrs/revision/coloc/nhgri/results/

trait=$1
type=$2
n=$3

# Get coloc files
./get_coloc_trait_nhgri.sh ${trait} ${type} ${n}
./summarize_coloc_nhgri.sh ${trait} > ${OUTDIR}/${trait}_coloc_results.tab
