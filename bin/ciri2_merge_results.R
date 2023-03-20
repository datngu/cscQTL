#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)
syntax='\nUsage:\t./ciri2_merge_results.R out_file libsize bsj_filter exp_prop in_file_1 in_file_2 in_file_3 ...\n\n'

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
    ciri2 = read.delim(file)
    ciri2 = ciri2[ciri2[,5] >= bsj_filter,]
    ciri2$bed_ID = paste(ciri2$chr, ciri2$circRNA_start, ciri2$circRNA_end, ciri2$gene_id, sep = "__")
    ID = unique( c(ID, ciri2$bed_ID))
}

all_circRNA = all_circRNA_count = data.frame(ID)

for(file in in_files){
    ciri2 = read.delim(file)
    #ciri2 = ciri2[ciri2[,5] >= bsj_filter,]
    s = basename(file)
    ciri2$bed_ID = paste(ciri2$chr, ciri2$circRNA_start, ciri2$circRNA_end, ciri2$gene_id, sep = "__")
    ciri2$TPM = ciri2[,5]/ sizes[s] * 10^6
    all_circRNA[,s] = ciri2$TPM[match(all_circRNA$ID, ciri2$bed_ID)]
    all_circRNA_count[,s] = ciri2[,5][match(all_circRNA$ID, ciri2$bed_ID)]
}

all_circRNA[is.na(all_circRNA)] = 0
all_circRNA_count[is.na(all_circRNA_count)] = 0

unexp_count <- apply(all_circRNA_count[,-1], 1, function(x) sum(x< bsj_filter))

all_circRNA = all_circRNA[unexp_count < length(in_files)*(1-exp_prop), ]

write.table(all_circRNA, file = out_file, row.names = F, col.names = T, quote = F, sep = "\t")
