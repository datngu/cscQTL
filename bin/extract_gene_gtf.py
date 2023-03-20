#!/usr/bin/env python

import sys
import gzip

in_fn = sys.argv[1]
out_fn = sys.argv[2]

print(in_fn)
print(out_fn)

#in_fn = "Salmo_salar-GCA_905237065.2-2021_07-genes.gtf.gz"
#out_fn = "gene.txt"

out_file = open(out_fn, 'w')

f = open(in_fn,'rt')
for l in f:
    if l[0] == "#":
        continue
    tem = l.split("\t")
    if tem[2] == "gene":
        tem2 = tem[8].split(";")[0]
        chr_id = tem[0]
        start = tem[3]
        end = tem[4]
        length = str( abs(int(tem[4]) - int(tem[3]) +1))
        strand = tem[6]
        gene_ID = tem2.split(" ")[1]
        res = chr_id + "\t" + start + "\t" + end + "\t" + gene_ID +  "\t"  + length + "\t" + strand + "\n"
        out_file.writelines(res)
f.close()
out_file.close()