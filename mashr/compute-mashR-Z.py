#!/usr/bin/env python3

import os
import pandas as pd
import sys

workdir = sys.argv[1]
prefix = sys.argv[2]

betas = pd.read_csv(os.path.join(workdir, "output-%s"%prefix, "posterior_betas.tsv"), sep="\t", index_col=0)
ses = pd.read_csv(os.path.join(workdir, "output-%s"%prefix, "posterior_beta_ses.tsv"), sep="\t", index_col=0)

z = betas/ses
z.to_csv(os.path.join(workdir, "output-%s"%prefix, "zscores.tsv"), sep="\t")
