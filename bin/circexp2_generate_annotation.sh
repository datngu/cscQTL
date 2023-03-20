#!/bin/bash
gtf=$1
out_name=$2

gtfToGenePred -genePredExt $gtf genePredExt.txt
awk '{print $12,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' genePredExt.txt > $out_name
rm genePredExt.txt