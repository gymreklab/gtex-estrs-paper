#!/usr/bin/env python

import sys

ensgene = sys.argv[1]
col = int(sys.argv[2])

gdict = {}
with open(ensgene, "r") as f:
    for line in f:
        ens, g = line.strip().split()
        gdict[g] = ens

for line in sys.stdin.readlines():
    items = line.strip().split()
    gene = gdict.get(items[col-1], "NA")
    sys.stdout.write("\t".join(map(str, items+[gene]))+"\n")
