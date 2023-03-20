#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)
syntax='\nUsage:\t./covariate_processing.R meta.csv exp_quatile_normalized.tsv genotype_PCs.tsv number_hidden_factors out_fn\n\n'

args = commandArgs(trailingOnly = TRUE)

if(length(args) < 5 ){
  cat("\nInvalid arguments, Program stop! \n")
  #cat("please see: https://github.com/PMBio/peer/wiki/Tutorial to understand peer input requirements")
  cat(syntax)
  quit()
}


# meta_path = "/home/datn/github/nf-circall-qtl/non_pipeline_scripts/meta.csv"
# exp_path = "/sigma4/projects/nf-circall-qtl/40t-cell-results/exp_quatile_normalized.tsv"
# genotype_PCs_path = "/sigma4/projects/nf-circall-qtl/40t-cell-results/enotype_4PCs.tsv"
# n = 1
# out_fn = "/sigma4/covariates.tsv"

meta_path = args[1]
exp_path = args[2]
genotype_PCs_path = args[3]
n = as.numeric(args[4])
out_fn = args[5]

cat("\nmeta input: ", exp_path, "\n")
cat("\nexpression input: ", exp_path, "\n")
cat("\ngenotype_PCs input: ", genotype_PCs_path, "\n")
cat("\nnumber_hidden_factors: ", n, "\n")
cat("\nout file name: ", out_fn, "\n")

if(!require(peer)){
    cat("failed to call peer package! please check your installation!")
    quit()
}
if(!require(data.table)){
    cat("failed to call data.table package! please check your installation!")
    quit()
}

meta = fread(meta_path)
exp = fread(exp_path)
genotype_PCs = fread(genotype_PCs_path, header = F)


## processing
meta = as.data.frame(meta)
row.names(meta) = meta[,2]
meta = meta[,-c(1,2)]

genotype_PCs = as.data.frame(genotype_PCs)
row.names(genotype_PCs) = genotype_PCs[,2]
genotype_PCs = genotype_PCs[,-c(1,2)]
colnames(genotype_PCs) = paste0("PC_", c(1: ncol(genotype_PCs)))

exp = as.data.frame(exp)
row.names(exp) = exp[,4]
exp = exp[,-c(1:4)]
expr = t(exp)
#dim(expr)
model = PEER()
PEER_setPhenoMean(model, as.matrix(expr))
PEER_setNk(model, n)
PEER_update(model)
peer_res = PEER_getX(model)
peer_res = as.data.frame(peer_res)
row.names(peer_res) = row.names(expr)
colnames(peer_res) = paste0("peer_", c(1:n))

all_covariates = cbind(meta, genotype_PCs)
all_covariates = cbind(all_covariates, peer_res)
id = row.names(all_covariates)
all_covariates = cbind(id, all_covariates)
all_covariates_final = as.data.frame(t(all_covariates))
write.table(all_covariates_final, file = out_fn, sep = "\t", col.names = F, row.names = T, quote = F)

cat("DONE!")