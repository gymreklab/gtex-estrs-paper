#!/bin/bash

OUTDIR=/storage/mgymrek/gtex-estrs/revision/coloc/results/

trait=$1
type=$2
n=$3

# Get coloc files
./get_coloc_trait.sh ${trait} ${type} ${n}
./summarize_coloc.sh ${trait} > ${OUTDIR}/${trait}_coloc_results.tab
