#!/usr/bin/env python3

"""
Example command:
CHROM=2
TISSUE=Lung
/storage/mgymrek/workspace/gtex-estrs-paper/str-regression/LinRegAssociationTest_v3.py \
  --chrom ${CHROM} \
  --expr /storage/szfeupe/Runs/650GTEx_estr/Analysis_by_Tissue/Review_Rerun/${TISSUE}/Corr_Expr.csv \
  --exprannot /storage/mgymrek/workspace/gtex-estrs-paper/str-regression/gencode_gene_annotations_hg19.csv \
  --strgt /storage/mgymrek/gtex-estrs/nonlinear/fromRichard/multiallel.csv.gz \
  --distfromgene 100000 \
  --out outfile.tab \
  --tmpdir /tmp \
  --nonlinear
"""

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

def GetKeepSamplesAlleles(locus_str, MINCOUNT):
    """
    locus_str is indexed by sample. has columns A1 and A2
    Remove samples with alleles seen < MINGENOTYPE times
    """
    all_alleles = list(locus_str["A1"])+list(locus_str["A2"])
    unique_alleles = set(all_alleles)
    acounts = dict([(al, all_alleles.count(al)) for al in unique_alleles])
    counts = locus_str.apply(lambda x: min([acounts[x["A1"]], acounts[x["A2"]]]), 1)
    return [item for item in locus_str.index if counts[item]>=MINCOUNT]

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

def ZNorm(vals, m=None, sd=None):
    if m is None: m = np.mean(vals)
    if sd is None: sd = math.sqrt(np.var(vals))
    if sd == 0: return None
    return [(item-m)/sd for item in vals]

def QuadraticRegression(X, Y, norm=False, norm_alleles=False, translate_alleles=None, minsamples=0):
    """
    input data has columns A1 and A2
    Perform regression Y~b1(x1+x2)+b2(x1**2+x2**2)
    Return b1, b1_se, b1_p, b2, b2_se, b2_p, overall p (F test), bic
    """
    if norm_alleles:
        alleles = list(X["A1"])+list(X["A2"])
        X1 = ZNorm(X["A1"], m=np.mean(alleles), sd=math.sqrt(np.var(alleles)))
        X2 = ZNorm(X["A2"], m=np.mean(alleles), sd=math.sqrt(np.var(alleles)))
    else:
        X1 = X["A1"]
        X2 = X["A2"]
    if translate_alleles is not None:
        X1 = [item+translate_alleles for item in X1]
        X2 = [item+translate_alleles for item in X2]
    if norm:
        Y = ZNorm(Y)
    if X1 is None or X2 is None or Y is None: return [None]*8
    if np.var(X1)==0 or np.var(X2)==0: return [None]*8
    if len(X1)<minsamples: return [None]*8
    X_lin = [X1[i]+X2[i] for i in range(len(X1))]
    X_quad = [X1[i]**2+X2[i]**2 for i in range(len(X2))]
    X = np.array([X_lin, X_quad]).transpose()
    X = sm.add_constant(X)
    mod_ols = sm.OLS(Y, X, missing="drop")
    res_ols = mod_ols.fit()
    b1, b2 = res_ols.params[1:]
    b1_se, b2_se = res_ols.bse[1:]
    b1_p, b2_p = res_ols.pvalues[1:]
    f_pval = res_ols.f_pvalue
    bic = res_ols.bic
    return b1, b1_se, b1_p, b2, b2_se, b2_p, f_pval, bic

def LinearRegression(X, Y, norm=False, minsamples=0, splitalleles=False):
    """
    Perform linear regression, return beta, beta_se, p, bic
    If splitalleles, input data has columns A1 and A2. Need to combine
    """
    if splitalleles:
        X = X["A1"]+X["A2"]
    if norm:
        X = ZNorm(X)
        Y = ZNorm(Y)
    if X is None or Y is None: return None, None, None, None
    if np.var(X)==0: return None, None, None, None
    if len(X) < minsamples: return None, None, None, None
    X = sm.add_constant(X)
    mod_ols = sm.OLS(Y, X, missing='drop')
    res_ols = mod_ols.fit()
    pval = res_ols.pvalues[1]
    slope = res_ols.params[1]
    err = res_ols.bse[1]
    return slope, err, pval, res_ols.bic

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
    parser.add_argument("--nonlinear", help="Perform non-linear association tests", action="store_true")

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
        if item not in expr.index: samples_to_remove.append(item)
    for item in samples_to_remove: str_samples.remove(item)
    expr = expr.loc[str_samples,:]
    PROGRESS("There are %s samples"%len(str_samples))
    
    f = open(OUTFILE, "w")
    header_items = ["gene", "chrom", "str.id", "str.start", "n.miss", "beta", "beta.se", "p.linear", "bic.linear"]
    if args.nonlinear: header_items.extend(["b1","b1_se","b1_p","b2","b2_se","b2_p","p_quad","bic_quad","best_model","best_p"])
    f.write("\t".join(header_items)+"\n")
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
            if args.nonlinear:
                # Separate the 2 alleles in A1 and A2
                locus_str = pd.DataFrame({"gt": locus_str})
                locus_str.index = list(cis_strs.index)
                samples_to_keep = [str_samples[k] for k in range(len(str_samples)) if "NA" not in locus_str.loc[str_samples[k],"gt"]]
                locus_str = locus_str.loc[samples_to_keep]
                locus_str["A1"] = locus_str["gt"].apply(lambda x: int(x.split(",")[0]))
                locus_str["A2"] = locus_str["gt"].apply(lambda x: int(x.split(",")[1]))
                # Filter outlier alleles (seen < 3 times)
                samples_to_keep = GetKeepSamplesAlleles(locus_str, MINGENOTYPE)
                locus_str = locus_str.loc[samples_to_keep]
            else:
                # Remove samples with NA genotypes, or with rare genotypes
                samples_to_keep = [str_samples[k] for k in range(len(str_samples)) if locus_str[str_samples[k]] != "None"]
                locus_str = locus_str[samples_to_keep].astype('float') #locus_str.loc[samples_to_keep,:]
                keepgts = GetKeepGTs([float(item) for item in locus_str], MINGT)
                samples_to_keep = [samples_to_keep[k] for k in range(len(samples_to_keep)) \
                                   if locus_str[samples_to_keep[k]] in keepgts]
                locus_str = locus_str[samples_to_keep]
                # Ignore locus if number of #genotype <3
                if len(set(locus_str)) < MINGENOTYPE:
                    continue

            # Remove filtered samples from expression data
            locus_y = y.loc[samples_to_keep,:]
            
            # Run regression using different models
            beta, beta_se, p_linear, bic_linear = LinearRegression(locus_str, locus_y["expr"].values, norm=NORM, minsamples=MINSAMPLES, splitalleles=args.nonlinear)

            if args.nonlinear:
                min_len = min([min(locus_str["A1"]), min(locus_str["A2"])])
                b1, b1_se, b1_p, b2, b2_se, b2_p, p_quad, bic_quad = QuadraticRegression(locus_str, locus_y["expr"].values, norm=NORM, \
                                                                                         translate_alleles=min_len, minsamples=MINSAMPLES)

            # Write output
            if beta is not None:
                items = [gene, CHROM, str_ids[j], starts[j], len(str_samples)-locus_str.shape[0], beta, beta_se, p_linear, bic_linear]
                if args.nonlinear:
                    if bic_quad is not None and bic_quad<bic_linear:
                        best_model = "QUAD"
                        best_p = p_quad
                    else:
                        best_model = "LINEAR"
                        best_p = p_linear
                    items.extend([b1, b1_se, b1_p, b2, b2_se, b2_p, p_quad, bic_quad, best_model, best_p])
                f.write("\t".join([str(item) for item in items])+"\n")
                if best_model == "QUAD" and p_quad<0.0001:
                    sys.stdout.write("\t".join([str(item) for item in items])+"\n") # TODO remove debug
    f.close()
