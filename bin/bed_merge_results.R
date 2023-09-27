#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)
syntax='\nUsage:\t./bed_merge_results.R out_file bsj_filter exp_prop in_file_1 in_file_2 in_file_3 ...\n\n'

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
# in_files = list.files( pattern = ".bed")
# collect sample IDs

ID = c()

for(file in in_files){
  mybed = fread(file)
  mybed = mybed[mybed$V6 >= bsj_filter,]
  mybed$bed_ID = paste(mybed$V1, mybed$V2, mybed$V3, mybed$V4, sep = "__")
  ID = unique( c(ID, mybed$bed_ID))
}


all_circRNA = all_circRNA_count = data.frame(ID)


for(file in in_files){
  mybed = fread(file)
  mybed$bed_ID = paste(mybed$V1, mybed$V2, mybed$V3, mybed$V4, sep = "__")

  s = basename(file)
  s = unlist(strsplit(s,".", fixed = TRUE))[1]
  mybed$TPM = mybed$V5
  mybed$bsj_count = mybed$V6
  #mybed$TPM = mybed$bsj_count
  all_circRNA[,s] = mybed$TPM[match(all_circRNA$ID, mybed$bed_ID)]
  all_circRNA_count[,s] = mybed$bsj_count[match(all_circRNA$ID, mybed$bed_ID)]

}


all_circRNA[is.na(all_circRNA)] = 0
all_circRNA_count[is.na(all_circRNA_count)] = 0

unexp_count <- apply(all_circRNA_count[,-1], 1, function(x) sum(x< bsj_filter))

all_circRNA = all_circRNA[unexp_count < length(in_files)*(1-exp_prop), ]

write.table(all_circRNA, file = out_file, row.names = F, col.names = T, quote = F, sep = "\t")
