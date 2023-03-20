#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)
syntax='\nUsage:\t./circall_merge_results.R out_file libsize bsj_filter exp_prop in_file_1 in_file_2 in_file_3 ...\n\n'

args = commandArgs(trailingOnly = TRUE)

if(length(args) < 4 ){
  cat("\nInvalid arguments, Program stop! \n")
  cat(syntax)
  quit()
}

out_file = args[1]
libsize_file = args[2]
bsj_filter = as.numeric(args[3])
exp_prop = as.numeric(args[4])
in_files = args[-c(1,2,3,4)]

cat("\nout file: ", out_file, "\n")
cat("\nlibsize_file: ", libsize_file, "\n")
cat("\nbsj_filter: ", bsj_filter, "\n")
cat("\nexp_prop: ", exp_prop, "\n")
x = paste(in_files, collapse = " ,")
cat("\nin file list:\n", x, "\n")



# sample_names = read.table(sample_file)
# sample_names = sample_names$V1

# collect sample IDs

libsize = read.delim(libsize_file, header = F)
sizes = libsize$V2
names(sizes) = libsize$V1

ID = c()
for(file in in_files){
    circall = read.delim(file)
    circall = circall[circall$junction_fragment_count >= bsj_filter,]
    circall$bed_ID = paste(circall$chr, circall$start, circall$end, circall$geneID, sep = "__")
    ID = unique( c(ID, circall$bed_ID))
}

all_circRNA = all_circRNA_count = data.frame(ID)

for(file in in_files){
  circall = read.delim(file)
  #circall = circall[circall$junction_fragment_count >= bsj_filter,]
  s = basename(file)
  circall$bed_ID = paste(circall$chr, circall$start, circall$end, circall$geneID, sep = "__")
  circall$TPM = circall$junction_fragment_count/ sizes[s] * 10^6
  all_circRNA[,s] = circall$TPM[match(all_circRNA$ID, circall$bed_ID)]
  all_circRNA_count[,s] = circall$junction_fragment_count[match(all_circRNA$ID, circall$bed_ID)]
}


all_circRNA[is.na(all_circRNA)] = 0
all_circRNA_count[is.na(all_circRNA_count)] = 0

unexp_count <- apply(all_circRNA_count[,-1], 1, function(x) sum(x< bsj_filter))

all_circRNA = all_circRNA[unexp_count < length(in_files)*(1-exp_prop), ]

write.table(all_circRNA, file = out_file, row.names = F, col.names = T, quote = F, sep = "\t")



