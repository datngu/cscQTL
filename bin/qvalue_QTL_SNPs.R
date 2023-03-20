#!/usr/bin/env Rscript

# Author: Dat T Nguyen <ndat@utexas.edu>
# Date: 24 Oct 2022

syntax="
Usage: this_script.R eqtl_normial=path/to/eqtl_normial.txt eqtl_permute=path/to/eqtl_permute.txt out_fn=path/to/out_file_name.txt
Input parameters are:
    eqtl_normial=path/to/eqtl_normial.txt
    eqtl_permute=path/to/eqtl_permute.txt
    out_fn=path/to/out_file_name.txt
"
#gwas_sumstat=path/to/sumstat.txt - should work with summary statistics downloaded from gwas catalog.
args = commandArgs(trailingOnly = TRUE)

options(stringsAsFactors = FALSE)
require(data.table)
require(qvalue)


out_fn = eqtl_normial = eqtl_permute = NA

if(length(args) < 3 ){
  cat("\nNot valid arguments, Program stop! \n")
  cat(syntax, useSource = TRUE)
  quit()
}

for (i in 1:length(args)){
  res=unlist(strsplit(args[i],"="))

  if (res[1] == "eqtl_normial") eqtl_normial = res[2]
  if (res[1] == "eqtl_permute") eqtl_permute = res[2]
  if (res[1] == "out_fn") out_fn = res[2]   
}

# setwd("/sigma4/projects/rnaQTL_hg38")
# eqtl_normial = "real_data_results/ciri2_qtl_mapping_nominal/nominal_pass_all_chrom.tsv.gz"
# eqtl_permute = "real_data_results/qtl_mapping_apply_qvalue/ciri2.tsv"
# out_fn = "test.txt"



do_qvalue <- function(gwas_sumstat, case_prop, qtl_all, qtl_hit){
    gwas_all = fread(gwas_sumstat)
    gwas_all = gwas_all[!is.na(gwas_all$p_value),]
    qtl_hit$PP.H4.abf=0
    
    res = list()
    for(i in 1:nrow(qtl_hit)){
        eqtl = qtl_all[qtl_all$V1 %in% qtl_hit$V1[i],]
        #eqtl$rs_id = eqtl$V2
        #eqtl$rs_id = eqtl$V8
        gwas = gwas_all[gwas_all$variant_id %in% eqtl$rs_id,]
        gwas$rs_id = gwas$variant_id
        gwas$p_value = as.numeric(gwas$p_value)
        input <- merge(eqtl, gwas, by="rs_id", all=FALSE, suffixes=c("_eqtl","gwas"))
        if(nrow(input) > 0){
        dataset1=list(pvalues=input$p_value, type="cc", s = case_prop, N=nrow(gwas))
        dataset2=list(pvalues=input$V12, type="quant", N=nrow(eqtl))
        result <- coloc.abf( dataset1, dataset2, MAF=input$effect_allele_frequency)
        #qtl_hit$PP.H3.abf[i] = result$summary[5]
        input_res = input[,c(2,1,13,18)]
        names(input_res) = c("ID", "rsID", "eqtl_pvalue", "gwas_pvalue" )
        n = as.character(input_res$ID[1])
        res[[n]] = input_res
        qtl_hit$PP.H4.abf[i] = result$summary[6]
        }
    }
    
    #out_name = paste0("coloc_result_", basename(gwas_sumstat))
    #fwrite(qtl_hit, file = out_name, sep = "\t")
    
    out_name2 = paste0("coloc_result_", basename(gwas_sumstat), ".Rdata")
    save(coloc_data = res, coloc_res = qtl_hit, file = out_name2)
}


qtl_all = fread(eqtl_normial, sep = " ")
qtl_all$rs_id = qtl_all$V8
qtl_all = qtl_all[!is.na(qtl_all$V12),]
qtl_hit = fread(eqtl_permute)

res = list()
for(i in 1:nrow(qtl_hit)){
    eqtl = qtl_all[qtl_all$V1 %in% qtl_hit$V1[i],]
    eqtl$qvalue = qvalue(eqtl$V12)$qvalues
    res[[i]] = eqtl
}

eqtl = qtl_all[qtl_all$V1 %in% qtl_hit$V1,]

qvalue = qvalue(eqtl$V12)

x =qvalue(eqtl$V12)$qvalue