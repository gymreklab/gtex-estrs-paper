#!/bin/R

#prefix='/storage/szfeupe/Runs/650GTEx_estr/Analysis_by_Tissue/'
#suffix='/Master.table'
prefix='temp/'
suffix='.table'
for (file in c('Artery-Aorta', 'Artery-Tibial', 'Lung')) {
	df = read.table(paste(prefix, file, suffix, sep=''), header=TRUE)
	print(file)
	print(sum(df$llqvalue < 0.05))
}
