#!/usr/bin/env python

import sys
import glob
import os
import argparse
import gzip
import pandas as pd

parser = argparse.ArgumentParser(description = "Meging circular dectected by Circall, CIRI2, CIRCexplorer2.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--gtf", help = "gtf file input")
parser.add_argument("--out", help = "output file name", default = "circRNA_consensus.txt")
parser.add_argument("--circall", nargs='+', help = "list of Circall results")
parser.add_argument("--ciri2", nargs='+', help = "list of CIRI2 results")
parser.add_argument("--circexp", nargs='+', help = "list of CIRCexplorer2 results")
parser.add_argument('--cutoff', type=int, default = 2, help = 'BSJ read count cutoff for filtering')
parser.add_argument('--consensus', type=int, default = 2, help = 'consensus cutoff: 2 means at least 2 methods called the targeted circRNAs')

args = parser.parse_args()

#print(args.gtf)

#print(args.circall)

#print(args.ciri2)

#print(args.circexp)

def is_gzipped(path):
    return path[-3:] == ".gz"


def get_gene_info(gtf):
    strand_dict = dict()
    if(is_gzipped(gtf)):
        f = gzip.open(gtf, 'rt')
    else:
        f = open(gtf, 'r')
    for l in f:
        if l[0] == "#":
            continue
        tem = l.split("\t")
        if tem[2] == "gene":
                tem2 = tem[8].split(";")[0]
                gene_ID = tem2.split(" ")[1].strip('\"')
                strand_dict[gene_ID] = tem[6]
    return strand_dict


# strand_dict = get_gene_info(args.gtf)

# for i in strand_dict:
#    print(i, "\t", strand_dict[i])


def read_circall(gtf, path_list, cutoff = 2):
    strand_dict = get_gene_info(gtf)
    res = pd.DataFrame()
    for path in path_list:
        # path = "/Users/datn/DATA_ANALYSES/circQTL_paper/nextflow_output/circall/EGAR00001193012"
        tem = pd.read_csv(path, sep = "\t", header=0)
        tem.start = tem.start
        tem.circID = tem.chr.astype(str) + "__" + tem.start.astype(str) + "__" + tem.end.astype(str) + "__" + tem.geneID
        tem = tem[tem.junction_fragment_count >= cutoff]
        tem["strand"] = pd.Series([strand_dict[x] for x in tem.geneID])
        tem = tem[["circID", "chr", "start", "end", "geneID", "strand"]]
        res = pd.concat([res,tem])
    res = res.drop_duplicates("circID")
    return res


def read_ciri2(path_list, cutoff = 2):
    res = pd.DataFrame()
    for path in path_list:
        # path = "/Users/datn/DATA_ANALYSES/circQTL_paper/nextflow_output/ciri2/EGAR00001193012"
        tem = pd.read_csv(path, sep = "\t", header=0)
        tem.circRNA_start = tem.circRNA_start
        tem.gene_id = tem.gene_id.str.split(",", expand=True)[0]
        tem.gene_id = tem.gene_id.fillna("Not_Available")
        tem["circID"] = tem.chr.astype(str) + "__" + tem.circRNA_start.astype(str) + "__" + tem.circRNA_end.astype(str) + "__" + tem.gene_id
        tem = tem[tem['#junction_reads'] >= cutoff]
        tem = tem[["circID", "chr", "circRNA_start", "circRNA_end", "gene_id", "strand"]]
        tem.columns = ["circID", "chr", "start", "end", "geneID", "strand"]
        res = pd.concat([res,tem])
    res = res.drop_duplicates("circID")
    return res

  

def read_ce2(path_list, cutoff = 2):
    res = pd.DataFrame()
    for path in path_list:
        # path = "/Users/datn/DATA_ANALYSES/circQTL_paper/nextflow_output/circexp2/EGAR00001193012"
        tem = pd.read_csv(path, sep = "\t", header= None)
        tem["chr"] = tem[0]
        tem["start"] = tem[1]+1
        tem["end"] = tem[2]
        tem["geneID"] = tem[14]
        tem["strand"] = tem[5]
        tem["circID"] = tem.chr.astype(str) + "__" + tem.start.astype(str) + "__" + tem.end.astype(str) + "__" + tem.geneID
        tem = tem[tem[12] >= cutoff]
        tem = tem[["circID", "chr", "start", "end", "geneID", "strand"]]
        res = pd.concat([res,tem])
    res = res.drop_duplicates("circID")
    return res



# gtf="/Users/datn/GENOMES/human/Homo_sapiens.GRCh38.106.gtf.gz"
# fi = glob.glob("/Users/datn/DATA_ANALYSES/circQTL_paper/nextflow_output/circall/EGAR*")
# circall = read_circall(gtf,fi)
# fi = glob.glob("/Users/datn/DATA_ANALYSES/circQTL_paper/nextflow_output/ciri2/EGAR*")
# ciri2 = read_ciri2(fi)  
# fi = glob.glob("/Users/datn/DATA_ANALYSES/circQTL_paper/nextflow_output/circexp2/EGAR*")
# ce2 = read_ce2(fi)

circall = read_circall(args.gtf, args.circall, cutoff = args.cutoff)
ciri2 = read_ciri2(args.ciri2, cutoff = args.cutoff)
ce2 = read_ce2(args.circexp, cutoff = args.cutoff)

all_circ = pd.concat([circall, ciri2, ce2])
vc = all_circ.circID.value_counts()
vc = vc[vc >= args.consensus]
pick = all_circ.circID.isin(vc.index)
all_circ = all_circ[pick].drop_duplicates("circID")

all_circ.to_csv(args.out, sep = "\t", index = False)


# python3 test.py --gtf /Users/datn/GENOMES/human/Homo_sapiens.GRCh38.106.gtf.gz --circall /Users/datn/DATA_ANALYSES/circQTL_paper/nextflow_output/circall/EGAR* --ciri2 /Users/datn/DATA_ANALYSES/circQTL_paper/nextflow_output/ciri2/EGAR* --circexp /Users/datn/DATA_ANALYSES/circQTL_paper/nextflow_output/circexp2/EGAR*