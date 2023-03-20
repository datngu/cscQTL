#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)
syntax='\nUsage:\t./salmon_merge_results.R gtf_path out_file in_file_1 in_file_2 in_file_3 ...\n\n'

args = commandArgs(trailingOnly = TRUE)

if(length(args) < 3 ){
  cat("\nInvalid arguments, Program stop! \n")
  cat(syntax)
  quit()
}


# gtf_path = "annotation.gtf"
# out_file = "test_txi"
# in_files = list.files(path=".", pattern="EGA", all.files = FALSE, full.names = F)


gtf_path = args[1]
out_file = args[2]
in_files = args[-c(1,2)]


cat("\nGTF file: ", gtf_path, "\n")
cat("\nout file: ", out_file, "\n")
x = paste(in_files, collapse = " ,")
cat("\nin file list:\n", x, "\n")


if(!require(GenomicFeatures)){
    cat("failed to call GenomicFeatures! please check your installation!")
    quit()
}
if(!require(tximport)){
    cat("failed to call tximport! please check your installation!")
    quit()
}

txdb = makeTxDbFromGFF(gtf_path)
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

names(in_files) = basename(in_files)
txi_salmon <- tximport(in_files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
txi_abundance = as.data.frame(txi_salmon$abundance)
ID = as.character(row.names(txi_abundance))
txi_abundance = cbind(ID, txi_abundance)
txi_abundance$ID = as.character(txi_abundance$ID)
row.names(txi_abundance) = NULL
write.table(txi_abundance, file = out_file, row.names = F, col.names = T, quote = F, sep = "\t")


