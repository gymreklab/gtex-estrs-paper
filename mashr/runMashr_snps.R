#!/usr/bin/env Rscript

chunksize = 10000

library('purrr')
library(mashr)
library(foreach)
library(doMC)
library(sqldf)

args = commandArgs(trailingOnly=TRUE)

if(length(args)==0){
    print("No arguments supplied.")
    stop("specify chrom")
}else{
    for(i in 1:length(args)){
         eval(parse(text=args[[i]]))
    }
}

workdir = '/storage/mgymrek/gtex-estrs/revision/mashr/'

indir = paste(workdir, '/input-snps-bychrom/chr', chrom, sep='')
outdir = paste(workdir, '/output-snps/chr', chrom, sep='')
intermediate = paste(workdir, '/intermediate-snps-bychrom/chr', chrom, sep='')
modeldir = paste(workdir, '/intermediate-snps-bychrom/chr1', sep='')

loadData = function(indir, intermediate) {
    print('----Load the data into two dataframes----')
    
    #get files to load
    files = list.files(path=indir, pattern="*.table", full.names=TRUE)
    shortFileNames = purrr::map(files, function(name) basename(tools::file_path_sans_ext(name)))

    #load files individually
    dataFrameList = purrr::map(files, read.table, header=TRUE, sep="\t", colClasses=c('character', 'factor', 'character', 'numeric', 'numeric','numeric'))

    #some rows are duplicated for some reason, so remove them
    dataFrameList = purrr::map(dataFrameList, function(x) x[!duplicated(x), ])

    #append the file name to the columns named beta and beta.se
    append = function (name) {
        function (df, ext) {
            names(df)[names(df) == name] = paste(name, ext, sep='_')
            return(df)
        }
    }

    dataFrameList = purrr::map2(dataFrameList, shortFileNames, append('beta'))
    dataFrameList = purrr::map2(dataFrameList, shortFileNames, append('beta.se'))

    #merge the dataframes
    allData = purrr::reduce(dataFrameList, function(df1, df2) base::merge(df1, df2, by=c('gene', 'chrom', 'str.start'), all=TRUE))

    #move gene, str idetnifier info into rownames
    rownames(allData) = purrr::pmap(allData[c(1,2,3)], paste, sep='_')
    allData = allData[-(1:3)]

    #separate betas and beta.ses
    betas = allData[grepl('beta_', colnames(allData))]
    beta.ses = allData[grepl('beta.se', colnames(allData))]
                            
    colnames(betas) = purrr::map(colnames(betas), substring, 6)
    colnames(beta.ses) = purrr::map(colnames(beta.ses), substring, 9)
    
    #figure out which rows contain something significant
    sig_pval_thresh = 1e-5

    # get a vector with the same row names, but all falses
    sigRows = betas["WholeBlood"] 
    colnames(sigRows) = "Significant"
    sigRows[TRUE] = FALSE
    naRows = betas["WholeBlood"]
    colnames(naRows) = "NARows"
    naRows[TRUE] = FALSE
    for (colname in colnames(betas)) {
        beta.ses[colname] = sapply(beta.ses[colname], as.numeric) # why need sapply?
        zscores = abs(betas[colname]/beta.ses[colname])[,1]
        zscores[betas[colname] == 0] = 0
        oneSidedReversePVals = pnorm(zscores)
        twoSidedPVals = 2*(1 - oneSidedReversePVals)
        sigRows = sigRows | (twoSidedPVals < sig_pval_thresh )
        naRows = naRows | is.na(betas[colname])
    }
    sigRows = sigRows & (!naRows) # NA rows cannot be in sigRows                       

    # Set NAs to 0s but big stderr
    betas[is.na(betas)] = 0
    beta.ses[is.na(beta.ses)] = 10

    saveRDS(betas, paste(intermediate, '/betas.rds', sep=''))
    saveRDS(beta.ses, paste(intermediate, '/beta.ses.rds', sep=''))
    saveRDS(sigRows, paste(intermediate, '/sigRows.rds', sep='')) # save for reuse
    saveRDS(naRows, paste(intermediate, '/naRows.rds', sep='')) # save for reuse

    # Return output
    return(list(betas, beta.ses, sigRows, naRows))
}

runMashrChunks = function(betas, beta.ses, sample_corr, fittedG, outdir, chunksize) {
    print('---running mashr by chunk---')
    registerDoMC(20) # Number of cores to use
    numchunks = ceiling(nrow(betas)/chunksize)
    foreach (idx=1:numchunks) %dopar% {
        print(paste('starting chunk ', idx, sep=''))
        maxval = (idx*chunksize)
        if (maxval > nrow(betas)) {
            maxval = nrow(betas)
        }
        rows=((idx-1)*chunksize+1):maxval
        chunkBetas = betas[rows,]
        chunkBeta.ses = beta.ses[rows,]
        prep = list(Bhat=data.matrix(chunkBetas), Shat=data.matrix(chunkBeta.ses))
        chunkMashrData = mash_set_data(prep$Bhat, prep$Shat, V=sample_corr)
        runMashr(chunkMashrData, fittedG, outdir, idx)
    }
}

runMashr = function(mashrData, fittedG, outdir, name) {
	print(paste('----run mashr with outdir ', outdir, '----', sep=''))
	library(mashr)

	mashrOutput = mash(mashrData, g=fittedG, fixg=TRUE)
	saveRDS(mashrOutput, paste(outdir, '/mashrOutput.rds', sep=''))

	lfsr = get_lfsr(mashrOutput)
	write.table(lfsr, paste(outdir, '/posterior_lfsr_', name, '.tsv', sep=''), sep="\t", col.names=NA, quote=FALSE)
	posterior_betas = get_pm(mashrOutput)
	write.table(posterior_betas, paste(outdir, '/posterior_betas_', name, '.tsv', sep=''), sep="\t", col.names=NA, quote=FALSE)
	posterior_beta_ses = get_psd(mashrOutput)
	write.table(posterior_beta_ses, paste(outdir, '/posterior_beta_ses_', name, '.tsv', sep=''), sep="\t", col.names=NA, quote=FALSE)
	log10bf = get_log10bf(mashrOutput)
	rownames(log10bf) = rownames(posterior_beta_ses)
	write.table(log10bf, paste(outdir, '/posterior_log10bf_', name, '.tsv', sep=''), sep="\t", col.names=NA, quote=FALSE)
}

# Step 0: Load precomputed model
fittedG = readRDS(paste(modeldir, '/fittedG.rds', sep=''))
sample_corr = readRDS(paste(modeldir, '/sample_corr.rds', sep=''))

# Step 1: load data
l = loadData(indir, intermediate)
betas = l[[1]]
beta.ses = l[[2]]
sigRows = l[[3]]
naRows = l[[4]]
#betas = readRDS(paste(intermediate, '/betas.rds', sep=''))
#beta.ses = readRDS(paste(intermediate, '/beta.ses.rds', sep=''))
#sigRows = readRDS(paste(intermediate, '/sigRows.rds', sep=''))
#naRows = readRDS(paste(intermediate, '/naRows.rds', sep=''))

# Step 2: Prepare mashR data using precomputed model
prep = list(Bhat = data.matrix(betas), Shat = data.matrix(beta.ses))
mashrData = mash_set_data(prep$Bhat, prep$Shat, V=sample_corr)

# Step 3: Run mashR
runMashrChunks(betas, beta.ses, sample_corr, fittedG, outdir, chunksize) 

