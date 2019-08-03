#!/usr/bin/env python3

import os
import pandas as pd
import sys

workdir = sys.argv[1]
prefix = sys.argv[2]

ZTHRESH = 3

# Read in zvals
zvals = pd.read_csv(os.path.join(workdir, "output-%s"%prefix, "zscores.tsv"), sep="\t", index_col=0)

def maxabs(x):
    vals = [abs(item) for item in x]
    return max(vals)

# Summarize - per tissue
all_estrs_genelevel = set()
all_genes_genelevel = set()
all_estrs_locuslevel = set()
all_genes_locuslevel = set()

for tissue in zvals.columns:
    tvals = zvals[[tissue]].copy()
    tvals["gene"] = [item.split("_")[0] for item in tvals.index]
    tvals["ID"] = tvals.index
    # Locus level
    tvals_locsig = tvals[tvals[tissue].apply(abs)>=ZTHRESH]
    all_estrs_locuslevel = all_estrs_locuslevel.union(set(tvals_locsig["ID"]))
    all_genes_locuslevel = all_genes_locuslevel.union(set(tvals_locsig["gene"]))
    tvals_locsig[[tissue]].to_csv(os.path.join(workdir, "output-%s"%prefix, "sig-bytissue", "%s-locuslevel-estrs.tsv"%tissue), sep="\t", index=True)
    print("%s: %s locus-level eSTRs"%(tissue, tvals_locsig.shape[0]))
    # Gene level
    tdata = tvals.groupby("gene", as_index=False).agg({tissue: maxabs})
    tdata.columns = ["gene","max"]
    tvals = pd.merge(tvals, tdata, on=["gene"])
    tvals = tvals[tvals[tissue].apply(abs)==tvals["max"]]
    tvals_sig = tvals[tvals["max"]>=ZTHRESH]
    print("%s: %s gene-level eSTRs"%(tissue, tvals_sig.shape[0]))
    all_estrs_genelevel = all_estrs_genelevel.union(set(tvals_sig["ID"]))
    all_genes_genelevel = all_genes_genelevel.union(set(tvals_sig["gene"]))
    tvals_sig.index = tvals_sig["ID"]
    tvals_sig[[tissue]].to_csv(os.path.join(workdir, "output-%s"%prefix, "sig-bytissue", "%s-genelevel-estrs.tsv"%tissue), sep="\t", index=True)

# Summarize - overall
print("Sig genes (best per tissue for each gene): %s eSTRs, %s unique genes"%(len(all_estrs_genelevel), len(all_genes_genelevel)))
print("Sig genes (all sig): %s eSTRs, %s unique genes"%(len(all_estrs_locuslevel), len(all_genes_locuslevel)))

