#!/usr/bin/env python3

import argparse
import gzip
import math
import numpy as np
import os
from collections import Counter
import pandas as pd
import random
import shutil
import statsmodels.api as sm
import sys
import tabix

EXPRFILE = None
EXPRANNOTFILE = None
DISTFROMGENE = None
STRGTFILE = ""
OUTFILE = ""
TMPDIR = "/tmp/"
CHECKCHROM = False
PERMUTE_EXPR = False
DEBUG = False
MINSAMPLES = 0
MINGENOTYPE = 3

def PROGRESS(msg):
    sys.stderr.write("%s\n"%msg.strip())

def GetKeepGTs(gtvals, MINGT):
    gts = set(gtvals)
    keepgt = [item for item in gts if not np.isnan(item) and str(item) != "None" and gtvals.count(item)>=MINGT]
    return keepgt

def GetCisSTRs(gtfile, chrom, start, end):
    tb = tabix.open(gtfile)
    records = tb.query(chrom, start-1, end+1)
    data = []
    for record in records:
        data.append(record)
    if len(data) == 0: return None
    df = pd.DataFrame(data)
    gtcols = [item.decode('UTF-8') for item in (gzip.open(gtfile, "r").readline().strip().split())]
    df.columns = gtcols
    return df

def LoadSamples(gtfile):
    samples = [item.decode('UTF-8') for item in (gzip.open(gtfile, "r").readline().strip().split()[2:])]
    return samples

def ZNorm(vals):
    m = np.mean(vals)
    sd = math.sqrt(np.var(vals))
    if sd == 0: return None
    return [(item-m)/sd for item in vals]

def LinearRegression(X, Y, norm=False, minsamples=0):
    """
    Perform linear regression, return beta, beta_se, p
    """
    if norm:
        X = ZNorm(X)
        Y = ZNorm(Y)
        if X is None or Y is None: return None, None, None
        if np.var(X)==0: return None, None, None
        if len(X) <= minsamples: return None, None, None
    X = sm.add_constant(X)
    mod_ols = sm.OLS(Y, X, missing='drop')
    res_ols = mod_ols.fit()
    pval = res_ols.pvalues[1]
    slope = res_ols.params[1]
    err = res_ols.bse[1]
    return slope, err, pval

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="get datasets for LMM assocation tests with STRs and SNPs")
    parser.add_argument("--expr", help="Normalized expression residuals", type=str, required=True)
    parser.add_argument("--exprannot", help="Expression annotation file", type=str, required=True)
    parser.add_argument("--chrom", help="Restrict analysis to this chromosome", type=str, required=False)
    parser.add_argument("--distfromgene", help="Look at STRs/SNPs within this distance of gene boundaries", type=int, required=True)
    parser.add_argument("--strgt", help="File with genotypes. Must be tabix indexed", type=str, required=True)
    parser.add_argument("--out", help="Write data files to this file", type=str, required=True)
    parser.add_argument("--scriptdir", help="Directory containing scripts", type=str, required=False)
    parser.add_argument("--checkchrom", help="Only load exons for relevant chromosome", action="store_true")
    parser.add_argument("--tmpdir", help="Tmp directory", type=str)
    parser.add_argument("--permute", help="Permute expression values", action="store_true")
    parser.add_argument("--norm", help="Normalize genotypes before doing association", action="store_true")
    parser.add_argument("--min-samples", help="Require data for this many samples", type=int, default=0)
    parser.add_argument("--min-genotypes", help="Require this many genotypes per test", type=int, default=3)
    parser.add_argument("--mingt", help="Remove genotypes with fewer than this many samples", type=int, default=1)
    parser.add_argument("--debug", help="Print debug info", action="store_true")
    args = parser.parse_args()
    EXPRFILE = args.expr
    EXPRANNOTFILE = args.exprannot
    CHROM = args.chrom
    if "chr" not in str(CHROM): CHROM="chr%s"%CHROM
    DISTFROMGENE = args.distfromgene
    STRGTFILE = args.strgt
    OUTFILE = args.out
    NORM = args.norm
    MINSAMPLES = args.min_samples
    MINGENOTYPE = args.min_genotypes
    MINGT = args.mingt
    if args.scriptdir is not None: SCRIPTDIR = args.scriptdir
    if args.checkchrom: CHECKCHROM = True
    if args.tmpdir is not None: TMPDIR = args.tmpdir
    if args.permute: PERMUTE_EXPR = True
    if args.debug: DEBUG = True

    # Load expression values
    PROGRESS("Load expression")
    if CHECKCHROM:
        x = list(pd.read_csv(EXPRFILE, nrows=1).columns.values)
        x = [item for item in x if item == "Unnamed: 0" or CHROM in item]
        expr = pd.read_csv(EXPRFILE, usecols=x) # reading in all the columns takes way too much memory,
                                                # pull out only ones we need
    else:
        expr = pd.read_csv(EXPRFILE)
    if "Unnamed: 0" in expr.columns:
        expr.index = expr["Unnamed: 0"].values
        expr = expr.drop("Unnamed: 0", 1)

    # Load expression annotation
    PROGRESS("Load annotation")
    expr_annot = pd.read_csv(EXPRANNOTFILE)
    expr_annot.index = expr_annot["probe.id"].values
    expr_annot = expr_annot.loc[[item for item in expr.columns if item in expr_annot.index],:]
    expr_annot = expr_annot[expr_annot["gene.chr"] == CHROM]

    # Restrict to STR samples
    str_samples = LoadSamples(STRGTFILE)
    samples_to_remove = []
    for item in str_samples:
        if item not in expr.index: samples_to_remove.append(item) #str_samples.remove(item)
    for item in samples_to_remove: str_samples.remove(item)
    expr = expr.loc[str_samples,:]
    PROGRESS("There are %s samples"%len(str_samples))
    
    f = open(OUTFILE, "w")
    f.write("\t".join(["gene", "chrom", "str.id", "str.start", "n.miss", "beta", "beta.se", "lambda.remel","p.wald"])+"\n")
    # For each gene:
    # Pull out STRs within distance of gene ends
    # For each STR - gene pair, get expr, str for samples with data and do linreg
    PROGRESS("Expression annotation size %s "%str(expr_annot.shape))
    for i in range(expr_annot.shape[0]):
        gene = expr_annot.index.values[i]
        PROGRESS(" Getting data for %s"%gene)
        start = expr_annot["gene.start"].values[i]
        end = expr_annot["gene.stop"].values[i]
        cis_strs = GetCisSTRs(STRGTFILE, CHROM, start-DISTFROMGENE, end+DISTFROMGENE)
        if cis_strs is None: continue
        PROGRESS("%s STRs tested \n"%str(cis_strs.shape[0]))

        # Preprocess all at once
        str_ids = ["STR_%s"%item for item in cis_strs["start"].values]
        starts = list(cis_strs["start"].values)
        cis_strs = cis_strs[str_samples].transpose()
        cis_strs.index = str_samples
        cis_strs.columns = str_ids

        if PERMUTE_EXPR:
            expr[gene] = random.sample(list(expr[gene].values), expr.shape[0])
        y = pd.DataFrame({"expr":list(expr.loc[:, gene])})
        y.index = str_samples

        for j in range(len(str_ids)):            
            # Get STR data
            locus_str = cis_strs.ix[:, str_ids[j]]

            # Filter
            samples_to_keep = [str_samples[k] for k in range(len(str_samples)) if locus_str[str_samples[k]] != "None"]
            locus_str = locus_str[samples_to_keep].astype('float') #locus_str.loc[samples_to_keep,:]
            locus_y = y.loc[samples_to_keep,:]
            
            
            # Remove samples with rare genotypes
            keepgts = GetKeepGTs([float(item) for item in locus_str], MINGT)
            samples_to_keep = [samples_to_keep[k] for k in range(len(samples_to_keep)) \
                               if locus_str[samples_to_keep[k]] in keepgts]
            locus_str = locus_str[samples_to_keep]
            locus_y = y.loc[samples_to_keep,:]
            
            #ignore locus if number of #genotype <3
            if len(set(locus_str)) < MINGENOTYPE:
                continue
          
            # Run regression
            beta, beta_se, p = LinearRegression(locus_str, locus_y["expr"].values, norm=NORM, minsamples=MINSAMPLES)
            # Write output
            if beta is not None:
                f.write("\t".join(map(str, [gene, CHROM, str_ids[j], starts[j], len(str_samples)-locus_str.shape[0], beta, beta_se, -1, p]))+"\n")
#        break
    f.close()
