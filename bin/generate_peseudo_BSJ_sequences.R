#!/usr/bin/env Rscript

# Author: Dat T Nguyen <ndat@utexas.edu>
# Date: 24 Oct 2022

syntax="
Usage: this_script.R genome=genome.fa circRNA=circRNA_info.txt out=out_file_name.fa
Input parameters are:
    genome=path/to/genome.fa
    circRNA=path/to/circRNA_info.txt
    out=path/to/out_file_name.fa
"

args = commandArgs(trailingOnly = TRUE)

options(stringsAsFactors = FALSE)
require(data.table)
require(Biostrings)



if(length(args) < 3 ){
  cat("\nNot valid arguments, Program stop! \n")
  cat(syntax, useSource = TRUE)
  quit()
}

for (i in 1:length(args)){
  res=unlist(strsplit(args[i],"="))

  if (res[1] == "genome") genome_path = res[2]
  if (res[1] == "circRNA") circrna_path = res[2]
  if (res[1] == "out") out = res[2]   
}



convertReverseComplement<-function(DNAseq){
  DNAarr=unlist(strsplit(DNAseq,""))
  #reverse
  DNAarr=rev(DNAarr)
  #complement
  Aid=which(DNAarr=="A")
  Tid=which(DNAarr=="T")
  Gid=which(DNAarr=="G")
  Cid=which(DNAarr=="C")
  DNAarr[Aid]="T"
  DNAarr[Tid]="A"
  DNAarr[Gid]="C"
  DNAarr[Cid]="G"
  #result
  DNAseqRc=paste(DNAarr,collapse = "")
  return(DNAseqRc) 
}

get_seq = function(genome, chr, s, e, strand){
  left_seq = substring(genome[[chr]], s, s+149)
  right_seq = substring(genome[[chr]], e-149, e)
  seq = paste0(right_seq, left_seq)
  if(strand == "-") seq = convertReverseComplement(seq)
  return(seq)
}


#genome_path = "/Users/datn/GENOMES/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
#circrna_path = "/Users/datn/github/rnaQTL_hg38/circRNA_consensus.txt"
#out="circRNA_consensus.fa"

genome = readDNAStringSet(genome_path)

chnames=sapply(names(genome), function(x) unlist(strsplit(x," "))[1])
chnames[which(chnames=="MT")]="M"
names(genome)=chnames

info = fread(circrna_path, header = T)
info$chr[which(info$chr=="MT")]="M"


## getting peseudo circRNA seqs
res = c()
for(i in 1:nrow(info)){
  chr = info$chr[i]
  s = info$start[i]
  e = info$end[i]
  strand = info$strand[i]
  t = get_seq(genome, chr, s, e, strand)
  ## avoid bugs when get_seq return empty
  if(length(t) == 1){
    names(t) = info$circID[i]
    res = c(res, t)
  }

}

seqs = DNAStringSet(res)
writeXStringSet(seqs, out)
  
  
  
#genome_path = "/Users/datn/GENOMES/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
#circrna_path = "/Users/datn/github/rnaQTL_hg38/circRNA_consensus.txt"
#out="circRNA_consensus.fa"
  
#./bin/generate_peseudo_BSJ_sequences.R genome=/Users/datn/GENOMES/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa circRNA=/Users/datn/github/rnaQTL_hg38/circRNA_consensus.txt out=circRNA_consensus.fa
  
  
    
  

