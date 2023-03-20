#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)
syntax='\nUsage:\t./quantile_norm_salmon.R salmon_exp meta_data gene_info exp_prop out_fn\n\n'
# output is qualtile normalized circular RNA expression, data formated as expression input for fastQTL tool.
# 1. rename samples to genotype sample names, according to meta data input.
# 2. quantile norm (log2+0.1) exp.

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
gene_info = args[3]
exp_prop = as.numeric(args[4])
out_fn = args[5]

cat("\nexpression data input: ", exp_path, "\n")
cat("\nmeta data input: ", meta_path, "\n")
cat("\ngene info: ", gene_info, "\n")
cat("\nexp prop: ", exp_prop, "\n")
cat("\nout file name: ", out_fn, "\n")

# ## functions
# quantile_norm <- function(df){
#   df_rank <- apply(df,2,rank,ties.method="min")
#   df_sorted <- data.frame(apply(df, 2, sort))
#   df_mean <- apply(df_sorted, 1, mean)
  
#   index_to_mean <- function(my_index, my_mean){
#     return(my_mean[my_index])
#   }
  
#   df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
#   rownames(df_final) <- rownames(df)
#   return(df_final)
# }

# exp_path = "/sigma4/projects/nf-rnaQTL/salmon/salmon_res_merged.tsv"
# meta_path = "/sigma4/dev/nf-rnaQTL/data/meta.csv"
# gene_info = "/sigma4/projects/nf-rnaQTL/bin/gene_info.tsv"
# exp_prop = 0.5
# out_fn = "test_out.tsv"

# loading
exp_raw = fread(exp_path)
exp_raw = as.data.frame(exp_raw)
row.names(exp_raw) = exp_raw$ID

meta = fread(meta_path)
meta = as.data.frame(meta)
meta = meta[meta$rna_id %in% colnames(exp_raw),]
rna_id = meta$rna_id
# rename rna_id to genotype id
exp_raw = exp_raw[,rna_id]
# NA to zero
exp_raw[is.na(exp_raw)] = 0
unexp_count <- apply(exp_raw, 1, function(x) sum(x<0.00000001))
exp_raw = exp_raw[unexp_count < length(rna_id)*(1-exp_prop), ]

## so, now we do z-scaling by genes
exp = as.matrix(exp_raw)
t_exp = t(exp)
t_exp_scaled = scale(t_exp, center = TRUE)
exp_scaled = t(t_exp_scaled)
# qunatile normalization by to have identical distribution by samples
# normalize.quantiles function is from package preprocessCore
exp = normalize.quantiles(exp_scaled)
exp = as.data.frame(exp)
names(exp) = meta$genotype_id
# exp_info = strsplit(exp_raw$ID,"__")
# exp_info = as.data.frame(do.call(rbind,exp_info))
# colnames(exp_info) = c("#Chr", "start",  "end")
# #exp_info$TargetID = exp_raw$ID
# exp_final = cbind(exp_info, exp)
# plot(density(exp[,12]))
# plot(density(exp[12,]))

## matching gene annotation
exp_info = fread(gene_info, header = F)
colnames(exp_info) = c("#Chr", "start", "end", "TargetID", "length", "strand")
exp_info = exp_info[, c("#Chr", "start", "end", "TargetID")]
exp_info = exp_info[ exp_info$TargetID %in% row.names(exp_raw),]
od = match(row.names(exp_raw), exp_info$TargetID)
exp_info = exp_info[od,]

exp_final = cbind(exp_info, exp)


fwrite(exp_final, out_fn, sep = "\t", row.names = F)
