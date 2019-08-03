#!/usr/bin/env python3
"""
Get overlap between FM-eSTRs and GWAS hits
"""

import argparse
import os
import pandas as pd
import pybedtools
import scipy.stats
import sys
import vcf

START_BUFFER = 5

def concat(x):
    return ",".join([str(item) for item in x])

def GetBiAllelicInt(x):
    if int(x) == 0: return 0
    else: return 1

def GetLD(str_record, snp_record, minmaf=0, mincount=3):
    sample_to_gts = {}
    all_str_alleles = set()
    allele_counts = {}
    allele_counts2 = {} # if using str2
    if None in str_record.ALT: allelelens = [0]
    else: allelelens = [0] + [len(item)-len(str_record.REF) for item in str_record.ALT]
    for sample in str_record:
        sample_to_gts[sample.sample] = {"STR": None, "SNP": None}
        if sample["GB"] is None: continue
        gb = sample["GB"].split("|")
        if "." in gb:
            continue
        else:
            alleles = [int(item) for item in gb]
        for a in alleles:
            all_str_alleles.add(a)
            allele_counts[a] = allele_counts.get(a, 0) + 1
        sample_to_gts[sample.sample]["STR"] = alleles
    for sample in snp_record:
        if sample.sample not in sample_to_gts.keys():
            continue
        if None not in sample.gt_alleles:
            sample_to_gts[sample.sample]["SNP"] = map(GetBiAllelicInt, sample.gt_alleles)
    str_data = []
    snp_data = []
    for sample in sample_to_gts:
        if sample_to_gts[sample]["STR"] is None or sample_to_gts[sample]["SNP"] is None: continue
        # Enforce allele counts
        if allele_counts[sample_to_gts[sample]["STR"][0]] < mincount: continue
        if allele_counts[sample_to_gts[sample]["STR"][1]] < mincount: continue
        # Get data
        str_data.append(sum(sample_to_gts[sample]["STR"]))
        snp_data.append(sum(sample_to_gts[sample]["SNP"]))
    return scipy.stats.pearsonr(str_data, snp_data)[0]

def CalcLD(snp_reader, str_reader, chrom, str_start, snp_start):
    str_record = None
    records = str_reader.fetch(chrom, str_start-START_BUFFER, str_start+START_BUFFER)
    for r in records:
        if r.INFO["START"] == str_start:
            str_record = r
            break
    if str_record is None:
        sys.stderr.write("ERROR: couldn't find STR record for %s:%s\n"%(chrom, str_start))
        return None
    records = snp_reader.fetch(chrom, snp_start-1, snp_start+1)
    snp_record = None
    for r in records:
        snp_record = r
        if snp_record.POS == snp_start: break
    if snp_record is None or abs(snp_record.POS-snp_start) >=2:
        sys.stderr.write("Could not find SNP locus %s:%s\n"%(chrom, snp_start))
        return None
    return GetLD(str_record, snp_record)
    
def GetR2(ld):
    if str(ld) != "None":
        return "%0.3f"%ld**2
    else: return "None"

def GetBestTissue(tissue_info):
    # Adipose-Subcutaneous_0.50_0.01_5.727238439277857e-16;Adipose-Visceral_0.42_0.01_2.5816312070383294e-08
    best = None
    best_pval = 1
    for tinfo in tissue_info.split(";"):
        tissue, beta, caviar, pval = tinfo.split("_")
        pval = float(pval)
        if pval < best_pval:
            best = tissue
            best_pval = pval
    return best

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--estrs", help="Table of eSTRs", type=str, required=True)
    parser.add_argument("--gwas", help="Table of GWAS hits", type=str, required=True)
    parser.add_argument("--window", help="Look for overlap in this window", type=int, default=50000)
    parser.add_argument("--caviar", help="CAVIAR threshold", type=float, default=0.3)
    parser.add_argument("--prefix", help="Prefix to name output", type=str, required=True)
    parser.add_argument("--outdir", help="Output directory", type=str, required=True)
    parser.add_argument("--str", help="VCF of STR genotypes", type=str, required=True)
    parser.add_argument("--snp", help="VCF of SNP genotypes", type=str, required=True)
    args = parser.parse_args()

    # Load eSTRs
    estrs = pd.read_csv(args.estrs, sep="\t")
    causal = estrs[estrs["score"]>=args.caviar]

    # Load GWAS
    gwas = pd.read_csv(args.gwas, sep="\t")

    # Get overlap within window
    gwas["window_start"] = (gwas["start"]-args.window).apply(lambda x: max([x, 0]))
    gwas["window_end"] = gwas["start"]+args.window
    if "chr" not in str(gwas["chrom"].values[0]): gwas["chrom"] = gwas["chrom"].apply(lambda x: "chr%s"%x)
    gwas["id"] = gwas.apply(lambda x: ":".join([str(x["chrom"]),str(x["start"]),x["rsid"],x["name"]]), 1)
    gwas_bed = pybedtools.BedTool.from_dataframe(gwas[["chrom","window_start","window_end","id"]])
    causal_bed = pybedtools.BedTool.from_dataframe(causal[["chrom","str.start","str.end"]])
    try:
        overlap = gwas_bed.intersect(causal_bed, wa=True, wb=True).to_dataframe()
    except pd.errors.EmptyDataError:
        sys.stderr.write("No overlaps found\n")
        sys.exit(1)
    overlap.columns = ["chrom","window_start","window_end","gwas_id","chrom_x","str.start","str.end"]
    overlap["gwas.rsid"] = overlap["gwas_id"].apply(lambda x: x.split(":")[2])
    overlap["gwas.pos"] = overlap["gwas_id"].apply(lambda x: int(x.split(":")[1]))
    overlap["gwas.trait"] = overlap["gwas_id"].apply(lambda x: x.split(":")[3])
    overlap = overlap[["chrom","str.start","gwas.rsid","gwas.pos","gwas.trait"]].drop_duplicates()

    # Candidates file for coloc
    cdata = pd.merge(overlap, causal[["chrom","str.start","gene.name","gene","tissue_info"]], on=["chrom","str.start"])
    cdata["rsid"] = cdata["gwas.rsid"]
    cdata["tissue"] = cdata["tissue_info"].apply(GetBestTissue)
    cdata = cdata.groupby(["gene","tissue","gene.name"], as_index=False).agg({"rsid": concat})
    cdata[["gene","tissue","gene.name","rsid"]].to_csv(os.path.join(args.outdir, args.prefix+"_candidates.tab"), sep="\t", index=False)

    # Put back eSTR data and collapse genes
    data = pd.merge(overlap, causal[["chrom","str.start","str.end", "score","tissue_info","str.motif.forward","str.motif.reverse","gene.name"]], on=["chrom","str.start"])
    data = data.groupby(["chrom","str.start", "str.end","gwas.rsid","gwas.pos","gwas.trait","str.motif.forward","str.motif.reverse"], as_index=False).agg({"gene.name": concat, "score": max})

    # Compute SNP-STR LD. Write dataset 3 as we go
    outf = open(os.path.join(args.outdir, args.prefix+"_SuppDataset.tsv"), "w")
    outf.write("\t".join(["chrom","str.start","str.end", "str.motif.forward","str.motif.reverse","score","gene.name","snp.str.ld","gwas.pos","gwas.rsid","trait"])+"\n")
    snp_reader = vcf.Reader(open(args.snp, "rb"))
    str_reader = vcf.Reader(open(args.str, "rb"))
    ldvals = []
    for i in range(data.shape[0]):
        ld = CalcLD(snp_reader, str_reader, data["chrom"].values[i].replace("chr",""), data["str.start"].values[i], data["gwas.pos"].values[i])
        items = [data["chrom"].values[i], data["str.start"].values[i], data["str.end"].values[i], data["str.motif.forward"].values[i], data["str.motif.reverse"].values[i], \
                 data["score"].values[i], data["gene.name"].values[i], GetR2(ld), \
                 data["gwas.pos"].values[i], data["gwas.rsid"].values[i], data["gwas.trait"].values[i]]
        line = "\t".join([str(item) for item in items])+"\n"
        sys.stdout.write(line)
        outf.write(line)
        ldvals.append(ld)
    data["snp.str.ld"] = ldvals
    outf.close()

    ##### Outputs ####
    # Separate table of things with exact overlap
    data["same"] = data.apply(lambda x: x["gwas.pos"]>=(x["str.start"]-2) and x["gwas.pos"]<=(x["str.end"]+2), 1)
    data[data["same"]][["chrom","str.start","str.end","str.motif.forward","str.motif.reverse","score","gene.name","gwas.pos","gwas.rsid","gwas.trait","snp.str.ld"]].to_csv(os.path.join(args.outdir, args.prefix+"_SuppDataset_overlap.tsv"), sep="\t", index=False)
    
    # How many total FM-eSTRs within 50kb?
    print("Total num STRs: %s"%(data[["chrom","str.start"]].drop_duplicates().shape[0]))
    # How many of these have r2>0.1
    print("Total num STRs r2>0.1: %s"%(data[data["snp.str.ld"]>0.1][["chrom","str.start"]].drop_duplicates().shape[0]))
    # How many of these have r2>0.8
    print("Total num STRs r2>0.8: %s"%(data[data["snp.str.ld"]>0.8][["chrom","str.start"]].drop_duplicates().shape[0]))
    # How many are exact same hit?
    print("Total num STRs exact hit: %s"%(data[data["same"]][["chrom","str.start"]].drop_duplicates().shape[0]))

if __name__ == "__main__":
    main()
