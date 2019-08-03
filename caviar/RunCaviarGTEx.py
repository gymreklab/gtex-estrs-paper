#!/usr/bin/env python3

# Example
# ./RunCaviarGTEx.py --zsnp /storage/mgymrek/gtex-estrs/revision/mashr/output-snps/zscores.tsv --zstr /storage/mgymrek/gtex-estrs/revision/mashr/output-strs/sig-bytissue/WholeBlood-estrs.tsv --tissue WholeBlood --samples /storage/mgymrek/gtex-estrs/revision/samples/WholeBlood.samples  --strgt /storage/mgymrek/gtex-estrs/revision/genotypes/GTExNormalizedSTRGenotypes.table.gz --snpgt /storage/mgymrek/gtex-estrs/revision//genotypes/GTExNormalizedSNPGenotypes_chr21.table.gz --out test.tab --genes ENSG00000160213.5 --tmpdir test/ --num-causal 2 --use-topn-snps 10

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import argparse
import gzip
import numpy as np
import math
import os
import pandas as pd
from subprocess import Popen, PIPE, DEVNULL
import sys
import tabix
import statsmodels.api as sm

def PROGRESS(msg, printit=True):
    if printit: # false for some messages when not in debug mode
        sys.stderr.write("%s\n"%msg.strip())

def LoadSamples(samplesfile):
    if samplesfile is None: return []
    return [item.strip() for item in open(samplesfile, "r").readlines()]

# Quickly load for small gene sets
def LoadReg(zfile, tissue, zthresh, genes, prefix="", tempdir="/tmp"):
    if len(genes) > 0:
        PROGRESS("Subsetting reg file to %s genes\n"%len(genes))
        # Write gene list to a file
        with open(os.path.join(tempdir, "genelist.txt"), "w") as f:
            for gene in genes: f.write(gene+"\n")
        # Grep for the gene or gene list from the linreg file
        newlinreg = os.path.join(tempdir, os.path.basename(zfile))
#        if os.path.exists(newlinreg):
#            PROGRESS("Warning: Intermediate linreg file %s exists. Overwriting\n"%os.path.join(tempdir, os.path.basename(zfile)))
        cmd1 = "head -n 1 %s > %s"%(zfile, newlinreg)
        cmd2 = "grep -f %s %s >> %s"%(os.path.join(tempdir, "genelist.txt"), \
                                      zfile, newlinreg)
        os.system(cmd1+";"+cmd2)
        # Open the new reduced regression file
        fname = newlinreg
    else:
        fname = zfile
    reg = pd.read_csv(fname, sep="\t", index_col=0)
    reg = reg[[tissue]]
    if abs(zthresh) > 0:
        reg = reg[abs(reg[tissue])>=abs(zthresh)]
    reg.columns = ["Z"]
    reg["absZ"] = reg["Z"].apply(abs)
    reg["str.start"] = [int(item.split("_")[-1]) for item in reg.index]
    reg["chrom"] = [item.split("_")[1] for item in reg.index]
    reg["gene"] = [item.split("_")[0] for item in reg.index]
    reg["ID"] = reg.apply(lambda x: prefix+x["chrom"]+":"+str(x["str.start"]), 1)
    return reg

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

def GetZ(xvals, yvals):
    X = np.array(xvals)
    Y = np.array(yvals)
    X = sm.add_constant(X)
    mod_ols = sm.OLS(Y, X, missing="drop")
    res_ols = mod_ols.fit()
    return res_ols.params[1]/res_ols.bse[1]

def RecomputeZ(str_genotypes, snp_genotypes, expr, mingt):
    str_zscores = []
    snp_zscores = []
    assert(len(str_genotypes)==1)
    # Get set of samples to keep
    unique_gts = set(str_genotypes[0])
    keepgt = [gt for gt in unique_gts if (not np.isnan(gt) and str_genotypes[0].count(gt)>=mingt)]
    def KeepGT(gt, keepgt):
        if np.isnan(gt): return False
        else: return gt in keepgt
    keep_samples = [i for i in range(len(str_genotypes[0])) if KeepGT(str_genotypes[0][i], keepgt)]
    # Recompute STR zscore
    expr_y = [expr[k] for k in keep_samples]
    str_x = [str_genotypes[0][k] for k in keep_samples]
    str_zscores.append(GetZ(str_x, expr_y))
    # Recompute SNP zscores
    for i in range(len(snp_genotypes)):
        snp_x = [snp_genotypes[i][k] for k in keep_samples]
        if len(set(snp_x))==1: snp_zscores.append(np.nan)
        else: snp_zscores.append(GetZ(snp_x, expr_y))
    keep_loci = [0] + [(i+1) for i in range(len(snp_zscores)) if not np.isnan(snp_zscores[i])]
    return str_zscores, snp_zscores, keep_samples, keep_loci

# Write LDFILE, ZFILE
def GenerateCAVIARFiles(gene, strreg, snpreg, strgt, snpgt, \
                        use_topn_strs, use_topn_snps, \
                        str_gt_ind, snp_gt_ind, \
                        expr, recompute_z, mingt, \
                        tmpdir):
    if not os.path.exists(os.path.join(tmpdir, gene)): os.mkdir(os.path.join(tmpdir, gene))
    strdata = strreg[strreg["gene"]==gene].sort_values("absZ", ascending=False).head(use_topn_strs).sort_values("str.start")
    snpdata = snpreg[snpreg["gene"]==gene].sort_values("absZ", ascending=False).head(use_topn_snps).sort_values("str.start")
    if strdata.shape[0] == 0 or snpdata.shape[0] == 0:
        return False
    # 1. Get LDFILE for only that set of variants. Reset strdata and snpdata
    str_genotypes, strdata = LoadGenotypes(strgt, str_gt_ind, strdata)
    snp_genotypes, snpdata = LoadGenotypes(snpgt, snp_gt_ind, snpdata)
    all_genotypes = pd.DataFrame(str_genotypes + snp_genotypes).transpose()
    # If we need to recompute Z scores, do it here and keep track of which samples we're using
    if recompute_z:
        str_zscores, snp_zscores, sample_ind, loc_ind = RecomputeZ(str_genotypes, snp_genotypes, \
                                                          list(expr.loc[gene]), mingt)
        all_genotypes = all_genotypes.iloc[sample_ind, loc_ind]
        strdata["Z"] = str_zscores
        snpdata["Z"] = snp_zscores
        snpdata = snpdata.iloc[[i-1 for i in loc_ind if i!=0],:] 
    # Z-normalize genotypes and compute LD
    for i in range(all_genotypes.shape[1]):
        all_genotypes.iloc[:,i] = ZNorm(all_genotypes.iloc[:,i])
    ldmatrix = (all_genotypes.corr())
    ldfile = os.path.join(tmpdir, gene, "LDFILE")
    ldmatrix.to_csv(ldfile, header=None, index=False, sep="\t")
    # 2. Get ZFILE
    zfile = os.path.join(tmpdir, gene, "ZFILE")
    strdata[["ID", "Z"]].to_csv(open(zfile, "w"), header=None, index=False, sep="\t")
    snpdata[["ID", "Z"]].to_csv(open(zfile, "a"), header=None, index=False, sep="\t")
    return True

def GetFloat(value):
    if value == "None": return np.nan
    else: return float(value)

def LoadGenotypes(gtfile, gtind, regdata):
    chrom = regdata["chrom"].values[0]
    start = min(regdata["str.start"])
    end = max(regdata["str.start"])
    positions = list(regdata["str.start"])
    loaded_positions = []
    tb = tabix.open(gtfile)
    records = tb.query(chrom, start-1, end+1)
    data = []
    for record in records:
        pos = int(record[1])
        if pos not in positions: continue
        loaded_positions.append(pos)
        data.append([GetFloat(record[i+2]) for i in gtind]) # first two cols are chrom, start
#    assert(len(positions)==len(loaded_positions))
#    assert([positions[i]==loaded_positions[i] for i in range(len(positions))])
    regdata = regdata[regdata["str.start"].apply(lambda x: x in loaded_positions)]
    assert([regdata["str.start"].values[i]==loaded_positions[i] for i in range(len(loaded_positions))])
    return data, regdata

def GetGenotypeIndices(strgtfile, snpgtfile, samples):
    str_samples = [item.decode('UTF-8') for item in (gzip.open(strgtfile, "r").readline().strip().split()[2:])]
    snp_samples = [item.decode('UTF-8') for item in gzip.open(snpgtfile, "r").readline().strip().split()[2:]]
    use_samples = list((set(str_samples).intersection(set(snp_samples))).intersection(samples))
    str_ind = [str_samples.index(item) for item in use_samples]
    snp_ind = [snp_samples.index(item) for item in use_samples]
    return str_ind, snp_ind, use_samples

# Run CAVIAR using LDFILE and ZFILE in tmp/
# Write output to tmp/
def RunCAVIAR(gene, tmpdir, numcausal):
    zfile = os.path.join(tmpdir, gene, "ZFILE")
    ldfile = os.path.join(tmpdir, gene, "LDFILE")
    outfile = os.path.join(tmpdir, gene, "CAVIAR")
    if os.path.exists(outfile+"_post"): os.remove(outfile+"_post")
    if os.path.exists(outfile+"_set"): os.remove(outfile+"_set")
    cmd = "CAVIAR -o %s -l %s -z %s -c %s"%(outfile, ldfile, zfile, numcausal)
    p = Popen(cmd, shell=True, stdout=DEVNULL, stderr=DEVNULL)
    output = p.communicate()[0]
    if p.returncode != 0:
        PROGRESS("CAVIAR on %s failed"%gene)
        return False
    return True

# Write output
# ["gene", "top_snp", "top_snp_score", "str", "str.score", "str.rank", "num.snps"]
def WriteOutput(outfile, gene, tmpdir):
    cav = pd.read_csv(os.path.join(tmpdir, gene, "CAVIAR_post"), sep="\t")
    cav.columns = ["variant", "x","posterior"]
    cav = cav.sort_values("posterior", ascending=False)
    cav["rank"] = [item+1 for item in range(cav.shape[0])]
    cav_snp = cav[cav["variant"].apply(lambda x: "SNP" in x)].sort_values("posterior", ascending=False)
    cav_str = cav[cav["variant"].apply(lambda x: "SNP" not in x)].sort_values("posterior", ascending=False)
    for i in range(cav_str.shape[0]):
        str_rank = cav_str["rank"].values[i]
        output = [gene, cav_snp["variant"].values[0], cav_snp["posterior"].values[0], \
                  cav_str["variant"].values[i], cav_str["posterior"].values[i], str_rank, \
                  cav_snp.shape[0]]
        output = [str(item) for item in output]
        outfile.write("\t".join(output)+"\n")
    outfile.flush()

# Write log of failed genes
def WriteLog(logfile, gene):
    logfile.write("Failed: %s\n"%gene)
    logfile.flush()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run CAVIAR on GTEx data")
    parser.add_argument("--zsnp", help="File with SNP zscores from mashR", type=str, required=True)
    parser.add_argument("--zstr", help="File with STR zscores from mashR", type=str, required=True)
    parser.add_argument("--tissue", help="Tissue to process", type=str, required=True)
    parser.add_argument("--samples", help="File with samples to process for this tissue", type=str, required=True)
    parser.add_argument("--strgt", help="File with noramlized STR genotypes", type=str, required=True)
    parser.add_argument("--snpgt", help="File with normalized SNP genotypes", type=str, required=True)
    parser.add_argument("--out", help="Write results to this file", type=str, required=True)
    parser.add_argument("--genes", help="Only process these genes", type=str, default="all")
    parser.add_argument("--genes-file", help="Only process genes in this file", type=str)
    parser.add_argument("--precomputed", help="Use precomputed input files for CAVIAR", action="store_true")
    parser.add_argument("--use-topn-strs", help="Use top n STRs (by p-value)", type=int, default=1000000)
    parser.add_argument("--use-topn-snps", help="Use top n SNPs (by p-value)", type=int, default=1000000)
    parser.add_argument("--zthresh", help="Only consider variants with at least this big of zscore", type=float, default=0)
    parser.add_argument("--mingt", help="Remove outlier STR samples. Requires use-topn-strs=1 and recompute-zscores=True", type=int, default=1)
    parser.add_argument("--recompute-z", help="Recompute Z scores using same set of samples for all variants. Requires use-topn-strs=1 and --expr file", action="store_true")
    parser.add_argument("--expr", help="Expression file", type=str, required=False)
    parser.add_argument("--num-causal", help="Number of causal variants to consider", type=int, default=3)
    parser.add_argument("--tmpdir", help="Use this directory for temporary files", type=str, default="/tmp")
    args = parser.parse_args()

    if (args.mingt > 1 and (args.use_topn_strs > 1 or not args.recompute_z)):
        PROGRESS("ERROR: if using --mingt you must have use-topn-strs=1 and recompute-z=True")
        sys.exit(1)
    if (args.recompute_z and args.use_topn_strs > 1):
        PROGRESS("ERROR: if using --recompute-z you must have use-topn-strs=1")
        sys.exit(1)

    # Get list of genes to process
    if args.genes != "all":
        genes = set(args.genes.split(","))
    elif args.genes_file is not None:
        genes = [item.strip() for item in open(args.genes_file, "r").readlines()]
    else:
        genes = []

    if not args.precomputed:
        PROGRESS("\nLoad strs regression")
        strreg = LoadReg(args.zstr, args.tissue, args.zthresh, genes, prefix="STR_", tempdir=args.tmpdir)

        # Load regression results
        PROGRESS("\nLoad snps regression")
        snpreg = LoadReg(args.zsnp, args.tissue, args.zthresh, genes, prefix="SNP_", tempdir=args.tmpdir)

        # Get list of samples to process
        samples = LoadSamples(args.samples)

        # Get sample indices for genotype data
        str_gt_ind, snp_gt_ind, samples = GetGenotypeIndices(args.strgt, args.snpgt, samples)

        if args.expr is not None:
            expr = pd.read_csv(args.expr).transpose()[samples]
        else: expr = None

        # Set genes list if not done already
        if len(genes) == 0: genes = set(strreg["gene"])

    # Prepare outputs
    outfile = open(args.out, "w")
    outfile.write("\t".join(["gene", "top_snp", "top_snp_score", "str", "str.score", "str.rank","num.snps"])+"\n")
    outfile.flush()
    logfile = open(args.out+".log", "w")
    logfile.flush()

    # For each gene:
    # 1. Get intermediate files
    # 2. Run CAVIAR
    # 3. Generate output
    for gene in genes:
        PROGRESS("Processing gene %s"%gene)
        if not args.precomputed: # Store in args.tmpdir/gene/ LDFILE, ZFILE
            if not GenerateCAVIARFiles(gene, strreg, snpreg, args.strgt, args.snpgt, \
                                       args.use_topn_strs, args.use_topn_snps, \
                                       str_gt_ind, snp_gt_ind, \
                                       expr, args.recompute_z, args.mingt, \
                                       args.tmpdir):
                WriteLog(logfile, gene)
                continue
        if not RunCAVIAR(gene, args.tmpdir, args.num_causal):
            WriteLog(logfile, gene)
            continue
        WriteOutput(outfile, gene, args.tmpdir)
    outfile.close()
    logfile.close()
