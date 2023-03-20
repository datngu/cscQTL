#!/usr/bin/env python

import sys
import os
import gzip

sample = sys.argv[1]
in_fn = sys.argv[2]
out_fn = sys.argv[3]

print(in_fn)
print(out_fn)


out_file = open(out_fn, 'w')
count = 0
with gzip.open(in_fn,'rt') as f:
    for l in f:
        count = count+1

count = count/4
out_line = sample + "\t" + str(count) + "\n"

out_file.writelines(out_line)
f.close()
out_file.close()