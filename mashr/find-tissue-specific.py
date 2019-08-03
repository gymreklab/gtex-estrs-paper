#!/usr/bin/env python3

import os
import pandas as pd
import sys

ZTHRESH = 4

try:
    zfile = sys.argv[1]
    outdir = sys.argv[2]
except:
    sys.stderr.write("Usage: ./find-tissue-specific.py <zscore>\n")
    sys.exit(1)

zscores = pd.read_csv(zfile, sep="\t", index_col=0)
tissues = list(zscores.columns)

def maxabs(vals):
    return max([abs(item) for item in vals])

for t in tissues:
    print(t)
    not_t = [item for item in tissues if item != t]
    tdata = zscores[[t]].copy()
    tdata["maxZ"] = zscores[not_t].apply(maxabs, 1)
    tdata["diff"] = tdata[t].apply(abs)-tdata["maxZ"]
    tdata["ID"] = tdata.index
    tdata = tdata.sort_values("diff", ascending=False)
    tdata[["ID",t,"maxZ","diff"]].to_csv(os.path.join(outdir, t+"_tissue_specific.tab"), sep="\t", index=False)

    print(tdata.head())
