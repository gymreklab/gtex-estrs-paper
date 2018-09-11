#!/usr/bin/env python

import sys
col=int(sys.argv[1])

for line in sys.stdin.readlines():
    items = line.strip().split()
    tdata = items[col-1].split(";")
    
    best_score = 0
    best_tissue = "NA"
    
    for item in tdata:
        tissue, beta, score = item.split("_")
        score = float(score)
        if score > best_score:
            best_score = score
            best_tissue = tissue

    items[col-1] = best_tissue
    sys.stdout.write("\t".join(map(str, items))+"\n")
