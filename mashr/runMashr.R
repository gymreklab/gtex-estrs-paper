#!/usr/bin/env Rscript

library('purrr')
library(mashr)
library(foreach)
library(doMC)
library(sqldf)

args = commandArgs(trailingOnly=TRUE)

if(length(args)==0){
    print("No arguments supplied.")
    runval = "test"
}else{
    for(i in 1:length(args)){
         eval(parse(text=args[[i]]))
    }
}

workdir = '/storage/mgymrek/gtex-estrs/revision/mashr/'

indir = paste(workdir, '/input-', runval, sep='')
outdir = paste(workdir, '/output-', runval, sep='')
intermediate = paste(workdir, '/intermediate-', runval, sep='')

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

prepMashr = function(betas, beta.ses, sigRows, naRows, intermediate) {
    print('----prep mashr----')
    library(mashr)
    #For overall coding pattern, look at this vignette:
    #https://stephenslab.github.io/mashr/articles/intro_mash.html

    #1)hand the data off to mashr
    prep = list(Bhat = data.matrix(betas), Shat = data.matrix(beta.ses))
    mashrData = mash_set_data(prep$Bhat[!naRows,], prep$Shat[!naRows,])
    sample_corr = estimate_null_correlation_simple(mashrData)
    mashrData = mash_update_data(mashrData, V=sample_corr)

    # Prepare other mashR matrices
    mashrDataStrong = mash_set_data(prep$Bhat[sigRows, ], prep$Shat[sigRows, ], V=sample_corr)
    mashrDataAll = mash_set_data(prep$Bhat, prep$Shat, V=sample_corr)

    #1.5)
    #account for correlation among samples
    saveRDS(sample_corr, paste(intermediate, '/sample_corr.rds', sep=''))

    #run pca then extreme deconvolution to learn effect patterns from the significant data
    if (runval == "test") {
      num_components = 2
    } else {
      num_components = 5
    }
    pca_cov_matrices = cov_pca(mashrDataStrong, num_components)
    ed_cov_matrices = cov_ed(mashrDataStrong, pca_cov_matrices) #?? what does this step do?
    canonical_cov_matrices = cov_canonical(mashrData)

    cov_matrices = c(canonical_cov_matrices, ed_cov_matrices)

    #3)fit the model and save its output
    mashrModel = mash(mashrData, cov_matrices, outputlevel = 1)
    fittedG = get_fitted_g(mashrModel)
    saveRDS(fittedG, paste(intermediate, '/fittedG.rds', sep=''))	

    return(list(mashrDataAll, fittedG, sample_corr))
}


runMashr = function(mashrData, fittedG, outdir) {
	print(paste('----run mashr with outdir ', outdir, '----', sep=''))
	library(mashr)

	mashrOutput = mash(mashrData, g=fittedG, fixg=TRUE)
	saveRDS(mashrOutput, paste(outdir, '/mashrOutput.rds', sep=''))

	lfsr = get_lfsr(mashrOutput)
	write.table(lfsr, paste(outdir, '/posterior_lfsr.tsv', sep=''), sep="\t", col.names=NA, quote=FALSE)
	posterior_betas = get_pm(mashrOutput)
	write.table(posterior_betas, paste(outdir, '/posterior_betas.tsv', sep=''), sep="\t", col.names=NA, quote=FALSE)
	posterior_beta_ses = get_psd(mashrOutput)
	write.table(posterior_beta_ses, paste(outdir, '/posterior_beta_ses.tsv', sep=''), sep="\t", col.names=NA, quote=FALSE)
	log10bf = get_log10bf(mashrOutput)
	rownames(log10bf) = rownames(posterior_beta_ses)
	write.table(log10bf, paste(outdir, '/posterior_log10bf.tsv', sep=''), sep="\t", col.names=NA, quote=FALSE)
}

runMashrChromByChrom = function(betas, beta.ses, sample_corr, fittedG, outdir) {
	print('---running mashr chrom by chrom---')
	#paralellize this process
	registerDoMC(23) #number of cores allotted

	chroms = paste("chr", c(1:22, 'X'), sep='')
	foreach (idx = 1:length(chroms)) %dopar% {
		chrom = chroms[idx]
		rows = grepl(paste(chrom, '_', sep=''), rownames(betas), fixed=TRUE)
		if (!any(rows)) {
			#in some tests, not all chromosomes will be present
			next
		}
		chromBetas = betas[rows, ]
		chromBeta.ses = beta.ses[rows, ]
		prep = list(Bhat = data.matrix(chromBetas), Shat = data.matrix(chromBeta.ses))
		chromMashrData = mash_set_data(prep$Bhat, prep$Shat, V=sample_corr)
		chromOutdir = paste(outdir, chrom, sep='/')
		dir.create(chromOutdir, showWarnings = FALSE)
		runMashr(chromMashrData, fittedG, chromOutdir)
	}
}

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

# Step 2: Prep mashR
l = prepMashr(betas, beta.ses, sigRows, naRows, intermediate)
mashrData = l[[1]]
fittedG = l[[2]]
sample_corr = l[[3]]

if (runval == "snps-bychrom/chr1") {
   stop("Done writing model. Not computing posteriors")
}

# Step 3: Run mashR
if (runval == "test") {
  runMashr(mashrData, fittedG, outdir)
} else {
  runMashrChromByChrom(betas, beta.ses, sample_corr, fittedG, outdir)
}
