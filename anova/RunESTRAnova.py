#!/usr/bin/env python3

# ./RunESTRAnova.py --sigsnps /storage/mgymrek/gtex-estrs/revision/mashr/output-snps/sig-bytissue/Adipose-Subcutaneous-estrs.tsv --sigstrs /storage/mgymrek/gtex-estrs/revision/mashr/output-strs/sig-bytissue/Adipose-Subcutaneous-estrs.tsv --samples /storage/mgymrek/gtex-estrs/revision/samples/Adipose-Subcutaneous.samples --strgt /storage/mgymrek/gtex-estrs/revision/genotypes/GTExNormalizedSTRGenotypes.table.gz --snpgt /storage/mgymrek/gtex-estrs/revision//genotypes/GTExNormalizedSNPGenotypes_chr21.table.gz --chrom chr21 --out test.tab --expr /storage/szfeupe/Runs/650GTEx_estr/Analysis_by_Tissue/Review_Rerun/Adipose-Subcutaneous/Corr_Expr.csv

import warnings
warnings.filterwarnings("ignore")

import argparse
import gzip
import math
import numpy as np
import os
import pandas as pd
import sys
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm
import tabix
import scipy.stats

"""
Analyze h2 gain in adding eSTRs vs. eSNPs using ANOVA
"""

def GetFloat(value):
    if value == "None": return np.nan
    else: return float(value)

def LoadEvars(sigfile, chrom):
    data = pd.read_csv(sigfile, sep="\t")
    data["chrom"] = data["ID"].apply(lambda x: x.split("_")[1])
    data["start"] = data["ID"].apply(lambda x: int(x.split("_")[2]))
    data["gene"] = data["ID"].apply(lambda x: x.split("_")[0])
    return data[["chrom","start","gene"]]
    
def LoadSamples(samplesfile):
    if samplesfile is None: return []
    return [item.strip() for item in open(samplesfile, "r").readlines()]

def GetGenotypeIndices(strgtfile, snpgtfile, samples):
    str_samples = [item.decode('UTF-8') for item in (gzip.open(strgtfile, "r").readline().strip().split()[2:])]
    snp_samples = [item.decode('UTF-8') for item in gzip.open(snpgtfile, "r").readline().strip().split()[2:]]
    use_samples = list((set(str_samples).intersection(set(snp_samples))).intersection(samples))
    str_ind = [str_samples.index(item) for item in use_samples]
    snp_ind = [snp_samples.index(item) for item in use_samples]
    return str_ind, snp_ind, use_samples

def PROGRESS(msg, printit=False):
    if printit: # false for some messages when not in debug mode
        sys.stderr.write("%s\n"%msg.strip())

def LoadGenotypes(gtfile, gtind, chrom, pos):
    start = pos-1
    end = pos+1
    tb = tabix.open(gtfile)
    records = tb.query(chrom, start-1, end+1)
    data = []
    for record in records:
        rpos = int(record[1])
        if rpos != pos: continue
        return [GetFloat(record[i+2]) for i in gtind]
    return None

def ZNorm(vals):
    vals = [float(x) for x in list(vals)]
    m = np.nanmean(vals)
    sd = math.sqrt(np.nanvar(vals))
    if sd == 0: return None
    if np.isnan(sd) or np.isnan(m): return None
    newvals = []
    for item in vals:
        if np.isnan(item):
            newvals.append(np.nan)
        else: newvals.append((item-m)*1.0/sd)
    return newvals

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyzing SNP+STR ANOVA")
    parser.add_argument("--sigsnps", help="File with best SNPs per gene", type=str, required=True)
    parser.add_argument("--sigstrs", help="File with best STRs per gene", type=str, required=True)
    parser.add_argument("--samples", help="File with samples to process for this tissue", type=str, required=True)
    parser.add_argument("--chrom", help="Process this chrom", type=str, required=True)
    parser.add_argument("--strgt", help="File with noramlized STR genotypes", type=str, required=True)
    parser.add_argument("--snpgt", help="File with normalized SNP genotypes", type=str, required=True)
    parser.add_argument("--out", help="Write results to this file", type=str, required=True)
    parser.add_argument("--mingt", help="Remove STR genotypes with fewer than this many samples", type=int, default=1)
    parser.add_argument("--expr", help="Normalized expression residuals", type=str, required=True)
    args = parser.parse_args()

    # Load sample info 
    samples = LoadSamples(args.samples)
    str_gt_ind, snp_gt_ind, samples = GetGenotypeIndices(args.strgt, args.snpgt, samples)

    # Load expression and annotation
    PROGRESS("Load expression")
    expr = pd.read_csv(args.expr).transpose()[samples]

    # Load eSNP/eSTR info
    PROGRESS("Load eSTR/eSNP info")
    esnps = LoadEvars(args.sigsnps, args.chrom)
    estrs = LoadEvars(args.sigstrs, args.chrom)
    genes = set(estrs[estrs["chrom"]==args.chrom]["gene"])

    outf = open(args.out, "w")
    outf.write("\t".join(["gene","STR","SNP","anova.pval"])+"\n")
    for gene in genes:
        PROGRESS("Processing %s..."%gene)
        strdata = estrs[estrs["gene"]==gene]
        str_pos = strdata["start"].values[0]
        snpdata = esnps[esnps["gene"]==gene]
        if snpdata.shape[0] == 0:
            PROGRESS("No eSNP for %s"%gene)
            outf.write("\t".join([gene, args.chrom+":"+str(str_pos), "NA","NA"])+"\n")
            continue
        snp_pos = snpdata["start"].values[0]
        str_genotypes = LoadGenotypes(args.strgt, str_gt_ind, args.chrom, str_pos)
        snp_genotypes = LoadGenotypes(args.snpgt, snp_gt_ind, args.chrom, snp_pos)
        if str_genotypes is None or snp_genotypes is None:
            PROGRESS("Error processing %s. Could not find STR or SNP genotypes"%gene)
            continue
        try:
            expr_vals = expr.loc[gene]
        except KeyError:
            PROGRESS("Error processing %s. No expression data present"%gene)
            continue            
        genedata = pd.DataFrame({"expr": expr_vals, "SNP": snp_genotypes, "STR": str_genotypes})
        genedata = genedata[~np.isnan(genedata["STR"]) & ~np.isnan(genedata["SNP"]) & ~np.isnan(genedata["expr"])]
        # Remove outlier STR genotypes
        gtcounts = genedata.groupby("STR", as_index=False).agg({"SNP": len})
        keepgt = set(gtcounts[gtcounts["SNP"]>=args.mingt]["STR"])
        genedata = genedata[genedata["STR"].apply(lambda x: x in keepgt)]
#        print("STR r: %s"%scipy.stats.pearsonr(genedata["STR"], genedata["expr"])[0])
#        print("SNP r: %s"%scipy.stats.pearsonr(genedata["SNP"], genedata["expr"])[0])
#        print(genedata.groupby("STR", as_index=False).agg({"SNP": len}))
        # Normalize
        genedata["STR"] = ZNorm(genedata["STR"])
        genedata["SNP"] = ZNorm(genedata["SNP"])
        genedata["expr"] = ZNorm(genedata["expr"])
        formula_snpstr = "expr ~ STR+SNP"
        formula_snp = "expr ~ SNP"
        try:
            lm_snpstr = ols(formula_snpstr, genedata).fit()
        except:
            PROGRESS("Error running snpstr model for gene: %s"%gene)
            continue
        try:
            lm_snp = ols(formula_snp, genedata).fit()
        except:
            PROGRESS("Error running SNP only model for gene: %s"%gene)
            continue            
        anova_results = anova_lm(lm_snp, lm_snpstr)
        pval = anova_results["Pr(>F)"].values[1]
        outitems = [gene, args.chrom+":"+str(str_pos), args.chrom+":"+str(snp_pos), pval]
        outf.write("\t".join([str(item) for item in outitems])+"\n")
