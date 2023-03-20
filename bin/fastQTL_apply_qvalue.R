#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)
syntax='\nUsage:\t./fastQTL_apply_qvalue.R out_file_prefix fdr in_file_1 in_file_2 in_file_3 ...\n\n'

args = commandArgs(trailingOnly = TRUE)

if(length(args) < 2 ){
  cat("\nInvalid arguments, Program stop! \n")
  cat(syntax)
  quit()
}

if(!require(R.utils)){
    cat("failed to call R.utils package! trying to install R.utils ...!")
    install.packages('R.utils', repos='http://cran.us.r-project.org', force = TRUE)
    require(R.utils)
}

if(!require(data.table)){
    cat("failed to call data.table package! trying to install data.table ...!")
    install.packages('data.table', repos='http://cran.us.r-project.org', force = TRUE)
    require(data.table)
}

if(!require(qvalue)){
    cat("failed to call qvalue package! trying to install qvalue ...!")
    if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos='http://cran.us.r-project.org')
    BiocManager::install("qvalue", force = TRUE)
    require(qvalue)
}




out_file = args[1]
fdr = as.numeric(args[2])
in_files = args[-c(1,2)]

cat("\nout file prefix: ", out_file, "\n")
cat("\nFDR:\n", fdr, "\n")
x = paste(in_files, collapse = " ,")
cat("\nin file list:\n", x, "\n")

# file = "13_peer_permute_pass_all_chrom.tsv.gz"
# in_files = c("10_peer_permute_pass_all_chrom.tsv.gz", "1_peer_permute_pass_all_chrom.tsv.gz", "11_peer_permute_pass_all_chrom.tsv.gz", "20_peer_permute_pass_all_chrom.tsv.gz", "12_peer_permute_pass_all_chrom.tsv.gz", "2_peer_permute_pass_all_chrom.tsv.gz", "13_peer_permute_pass_all_chrom.tsv.gz", "3_peer_permute_pass_all_chrom.tsv.gz", "14_peer_permute_pass_all_chrom.tsv.gz", "4_peer_permute_pass_all_chrom.tsv.gz", "15_peer_permute_pass_all_chrom.tsv.gz", "5_peer_permute_pass_all_chrom.tsv.gz", "16_peer_permute_pass_all_chrom.tsv.gz", "6_peer_permute_pass_all_chrom.tsv.gz", "17_peer_permute_pass_all_chrom.tsv.gz", "7_peer_permute_pass_all_chrom.tsv.gz", "18_peer_permute_pass_all_chrom.tsv.gz", "8_peer_permute_pass_all_chrom.tsv.gz", "19_peer_permute_pass_all_chrom.tsv.gz", "9_peer_permute_pass_all_chrom.tsv.gz")
# fdr = 0.05
# out_file = "circall"


res_all = list()
res_all_raw = list()

df = data.frame(matrix(ncol = 2, nrow = 0))

for(file in in_files){
  n_peer = gsub("_peer_permute_pass_all_chrom.tsv.gz", "", basename(file))
  d = fread(file)
  d=d[!is.na(d$V17),]
  d$qval=qvalue(d$V17)$qvalue
  res_all_raw[[n_peer]] = d
  x = d[d$qval< fdr,]
  res_all[[n_peer]] = x
  n_p = as.integer(n_peer)
  n_c = nrow(x)
  df = rbind(df, c(n_p, n_c))
}
colnames(df) = c("n_peer", "n_ecircQTL")

od = order(df$n_peer)
df = df[od,]

# out_file = "circall"
# Give the chart file a name
png(file = paste0(out_file, ".png"))
barplot(df$n_ecircQTL ,names.arg = as.character(df$n_peer), xlab="number of peer factors", ylab="number of ecircQTL",
        main = out_file)
dev.off()

# finding optimal number of peer factors
optimal_peer = 1
n = 0
for(i in 1:nrow(df)){
  if(df$n_ecircQTL[i] > n){
    n = df$n_ecircQTL[i]
    optimal_peer = i
  }
}
optimal_peer = as.character(optimal_peer)
ecircQTL = res_all[[optimal_peer]]
ecircQTL_raw = res_all_raw[[optimal_peer]]

# export output
fwrite(ecircQTL, file = paste0(out_file, ".tsv"), sep = "\t", col.names = F)
fwrite(ecircQTL_raw, file = paste0(out_file, "_raw.tsv"), sep = "\t", col.names = F)
fwrite(df, file = paste0(out_file, ".stat"), sep = "\t")
fwrite( as.list(optimal_peer), file = out_file)




