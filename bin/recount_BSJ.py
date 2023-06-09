#!/usr/bin/env python

import sys
import os
import argparse
import pandas as pd


parser = argparse.ArgumentParser(description = "Counting BSJ read output by Circall_BSJ", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--bsj", help = "Read mapping info generated by Circall_BSJ", required = True)
parser.add_argument("--info", help = "Fragment mapping info generated by Circall_wt - used for library size normalization", required = True)
parser.add_argument("--out", help = "Output file name", default = "circRNA_re_count.txt")
parser.add_argument('--read', default=True, action='store_true', help = "Counting BSJ reads; Mutually exclusive with --fragment")
parser.add_argument('--fragment', dest='read', default=False, action='store_false', help = "Counting BSJ fragments. Mutually exclusive with --read")
parser.add_argument('--hit', type=int, default = 70, help = 'Min hit length')
parser.add_argument('--anchor', type=int, default = 10, help = 'Min anchor length')

args = parser.parse_args()


#print(args.read)
bsj_fn = args.bsj
out_fn = args.out
info_fn = args.info
read = args.read
min_hit_len = args.hit
min_anchor_len = args.anchor

# min_hit_len = 50
# min_anchor_len = 10

# bsj_fn = "/sigma4/data/genome_ref_hg38/test/sample12.txt_tem_dir/BSJ_mapped_read.txt"
# out_fn = "/sigma4/data/genome_ref_hg38/test/sample12.txt_tem_dir/recount1.txt"
# info_fn = "/sigma4/data/genome_ref_hg38/test/sample12.txt_tem_dir/Circall_wt_fragmentInfo.txt"
# fragment = True

read_info = pd.read_csv(bsj_fn, header = None, sep = "\t")
# some processing
read_info.columns = ["Header", "Direction", "QueryPos", "HitLen", "HitPos", "Txname", "Reflen"]
read_info["dis5"] = read_info.Reflen/2 - read_info.HitPos
read_info["dis3"] = read_info.HitPos + read_info.HitLen - read_info.Reflen/2
read_info["read_mate"] = read_info.Header.str.split(" ", expand = True)[0]
read_info["read_flag"] = read_info.Header.str.split("___", expand = True)[1]
read_info["read_id_tx_id"] = read_info.read_mate + "__" + read_info.Txname
## do some filtering
read_info = read_info[read_info.HitLen >= min_hit_len]
read_info = read_info[read_info.dis5 >= min_anchor_len]
read_info = read_info[read_info.dis3 >= min_anchor_len]

## create count data

count_data = pd.DataFrame(read_info.Txname.unique())
count_data.columns = ["circRNA"]
count_data.index = count_data["circRNA"]
count_data["bsj_count"] = 0


if read == True:
    #read_info = read_info.drop_duplicates("read_id_tx_id")
    for i, tx in enumerate(read_info.Txname):
        count_data.loc[tx, "bsj_count"] +=1


if read == False:
    read_info = read_info.drop_duplicates("read_id_tx_id")
    for i, tx in enumerate(read_info.Txname):
        count_data.loc[tx, "bsj_count"] +=1


## libsize normalization
info = pd.read_csv(info_fn, header = 0, sep = "\t")
libsize = info.numObservedFragments[0]/2

count_data["TPM"] = count_data["bsj_count"]/libsize * 10**6

count_data.to_csv(out_fn, index=False, sep = "\t")




# ./recount_BSJ.py --bsj BSJ_mapped_read.txt --out recount1.txt --info Circall_wt_fragmentInfo.txt