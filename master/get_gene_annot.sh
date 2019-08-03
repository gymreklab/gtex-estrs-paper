#!/bin/bash

# Get gene, gene.name, gene.strand

GTF=/storage/resources/dbase/human/gene_annotations/gencode.v19.genes.v7.patched_contigs.gtf

echo "gene,gene.name,gene.strand,gene.type" | sed 's/,/\t/g'
cat $GTF | grep -w gene | \
    cut -f 7,9 | cut -d';' -f1,3,5 | \
    sed 's/gene_id //' | sed 's/ gene_type //' | sed 's/ gene_name //' | \
    sed 's/"//g' | sed 's/;/\t/g' | \
    awk '{print $2 "\t" $4 "\t" $1 "\t" $3}'
    
    
