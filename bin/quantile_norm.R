#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)
syntax='\nUsage:\t./quantile_norm.R circularRNA_exp meta_data out_fn\n\n'
# output is qualtile normalized circular RNA expression, data formated as expression input for fastQTL tool.
# 1. rename samples to genotype sample names, according to meta data input.
# 2. z-scaling and quantile norm 

args = commandArgs(trailingOnly = TRUE)

if(length(args) < 3 ){
  cat("\nInvalid arguments, Program stop! \n")
  cat(syntax)
  quit()
}

if(!require(data.table)){
    cat("failed to call data.table package! please check your installation!")
    quit()
}

if(!require(preprocessCore)){
    cat("failed to call preprocessCore package! please check your installation!")
    quit()
}

exp_path = args[1]
meta_path = args[2]
out_fn = args[3]

cat("\nexpression data input: ", exp_path, "\n")
cat("\nmeta data input: ", meta_path, "\n")
cat("\nout file name: ", out_fn, "\n")

# exp_path = "/sigma4/projects/nf-circall-qtl/40t-cell-results/circall_res_merged"
# meta_path = "/sigma4/dev/nf-rnaQTL/data/meta.csv"
#out_fn = "test_out.tsv"
# loading
exp_raw = fread(exp_path)
exp_raw = as.data.frame(exp_raw)
row.names(exp_raw) = exp_raw$ID
meta = fread(meta_path)
meta = as.data.frame(meta)
meta = meta[meta$rna_id %in% colnames(exp_raw),]
rna_id = meta$rna_id
exp_raw = exp_raw[,rna_id]
exp = as.matrix(exp_raw)

t_exp = t(exp)
t_exp_scaled = scale(t_exp, center = TRUE)
exp_scaled = t(t_exp_scaled)
# quantile normalization by to have identical distribution by samples
# normalize.quantiles function is from package preprocessCore
exp = normalize.quantiles(exp_scaled)
exp = as.data.frame(exp)
names(exp) = meta$genotype_id
exp_info = strsplit(row.names(exp_raw),"__")
exp_info = as.data.frame(do.call(rbind,exp_info))
exp_info = exp_info[,c(1:3)]
colnames(exp_info) = c("#Chr", "start", "end")
exp_info$TargetID = row.names(exp_raw)
exp_final = cbind(exp_info, exp)
fwrite(exp_final, out_fn, sep = "\t", row.names = F)