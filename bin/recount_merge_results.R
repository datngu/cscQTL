#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)
syntax='\nUsage:\t./circall_merge_results.R out_file bsj_filter exp_prop in_file_1 in_file_2 in_file_3 ...\n\n'

args = commandArgs(trailingOnly = TRUE)

if(length(args) < 4 ){
  cat("\nInvalid arguments, Program stop! \n")
  cat(syntax)
  quit()
}

out_file = args[1]
bsj_filter = as.numeric(args[2])
exp_prop = as.numeric(args[3])
in_files = args[-c(1,2,3)]

cat("\nout file: ", out_file, "\n")
cat("\nbsj_filter: ", bsj_filter, "\n")
cat("\nexp_prop: ", exp_prop, "\n")
x = paste(in_files, collapse = " ,")
cat("\nin file list:\n", x, "\n")


require(data.table)
# in_files = list.files( pattern = ".recount")
# collect sample IDs

ID = c()
for(file in in_files){
    circall = fread(file)
    circall = circall[circall$bsj_count >= bsj_filter,]
    info = strsplit(circall$circRNA, split = "__", fixed = T)
    circall$chr = sapply(info, FUN = function(x){unlist(x)[1]})
    circall$start = sapply(info, FUN = function(x){unlist(x)[2]})
    #circall$start = as.integer(circall$start)
    circall$end = sapply(info, FUN = function(x){unlist(x)[3]})
    circall$geneID = sapply(info, FUN = function(x){unlist(x)[4]})
    circall$bed_ID = paste(circall$chr, circall$start, circall$end, circall$geneID, sep = "__")
    ID = unique( c(ID, circall$bed_ID))
}

all_circRNA = all_circRNA_count = data.frame(ID)

for(file in in_files){
    circall = fread(file)
    #circall = circall[circall$bsj_count >= bsj_filter,]
    info = strsplit(circall$circRNA, split = "__", fixed = T)
    circall$chr = sapply(info, FUN = function(x){unlist(x)[1]})
    circall$start = sapply(info, FUN = function(x){unlist(x)[2]})
    #circall$start = as.integer(circall$start)
    circall$end = sapply(info, FUN = function(x){unlist(x)[3]})
    circall$geneID = sapply(info, FUN = function(x){unlist(x)[4]})
    circall$bed_ID = paste(circall$chr, circall$start, circall$end, circall$geneID, sep = "__")

    s = basename(file)
    s = unlist(strsplit(s,".", fixed = TRUE))[1]
    circall$TPM = circall$TPM
    #circall$TPM = circall$bsj_count
    all_circRNA[,s] = circall$TPM[match(all_circRNA$ID, circall$bed_ID)]
    all_circRNA_count[,s] = circall$bsj_count[match(all_circRNA$ID, circall$bed_ID)]

}

all_circRNA[is.na(all_circRNA)] = 0
all_circRNA_count[is.na(all_circRNA_count)] = 0

unexp_count <- apply(all_circRNA_count[,-1], 1, function(x) sum(x< bsj_filter))

all_circRNA = all_circRNA[unexp_count < length(in_files)*(1-exp_prop), ]

write.table(all_circRNA, file = out_file, row.names = F, col.names = T, quote = F, sep = "\t")
