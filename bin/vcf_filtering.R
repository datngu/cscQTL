#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)
syntax='\nUsage:\t./vcf_filtering.R in.vcf.gz meta_path maf_cutoff out.vcf.gz\n\n'
# output is filtered vcf.gz file.

args = commandArgs(trailingOnly = TRUE)

if(length(args) < 4 ){
  cat("\nInvalid arguments, Program stop! \n")
  cat(syntax)
  quit()
}

if(!require(data.table)){
    cat("failed to call data.table package! please check your installation!")
    quit()
}

vcf_in = args[1]
meta_path = args[2]
maf_cutoff = as.numeric(args[3])
vcf_out = args[4]

cat("\nvcf input: ", vcf_in, "\n")
cat("\nmeta data input: ", meta_path, "\n")
cat("\nmaf cutoff: ", maf_cutoff, "\n")
cat("\nvcf output: ", vcf_out, "\n")

# vcf_in = "/sigma4/data/40T-cell-blueprint/snp_data/tcel.vcf.gz"
# meta_path = "/home/datn/github/nf-circall-qtl/non_pipeline_scripts/meta.csv"
# maf_cutoff = 0.05
# vcf_out = "/sigma4/data/40T-cell-blueprint/snp_data/tcel_filtered.vcf.gz"
# require(data.table)


meta = fread(meta_path)
genotype_id = meta[,2]
fwrite(genotype_id, "genotype_id.txt", row.names = F, col.names = F)
major_cutoff = 1- maf_cutoff

cmd = paste0( "bcftools view ", vcf_in, " -m2 -M2 -v snps -Q ", major_cutoff, ":major -q ", maf_cutoff, ":minor -e \'ALT=\".\"\' -S genotype_id.txt -Oz -o ", vcf_out)
system(cmd)

cat("\nDONE!!\n")

#cmd = paste0("bcftool view ", vcf_in, " -s ", genotype_id, " ")

#bcftools view $in_ref  -S $ref_dir/${sample}/sampleID_by_line.txt -a -o $ref_dir/${sample}/${sample}_ref_panel.vcf.gz -Oz