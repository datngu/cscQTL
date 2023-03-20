#!/usr/bin/env Rscript

# Author: Dat T Nguyen <ndat@utexas.edu>
# Date: 15 Sep 2022

syntax="
Usage: coloc.R meta=path/to/meta.csv eqtl_normial=path/to/eqtl_normial.txt eqtl_permute=path/to/eqtl_permute.txt
Input parameters are:
    meta=path/to/meta.csv - meta file with 2 columns, the first is path to gwas sumstat, the second is proportion of case samples in the gwas.
    eqtl_normial=path/to/eqtl_normial.txt
    eqtl_permute=path/to/eqtl_permute.txt
    
"
#gwas_sumstat=path/to/sumstat.txt - should work with summary statistics downloaded from gwas catalog.
args = commandArgs(trailingOnly = TRUE)

options(stringsAsFactors = FALSE)
require(data.table)
require(coloc)

meta = eqtl_normial = eqtl_permute = gwas_sumstat = case_prop = NA

if(length(args) < 3 ){
  cat("\nNot valid arguments, Program stop! \n")
  cat(syntax, useSource = TRUE)
  quit()
}

for (i in 1:length(args)){
  res=unlist(strsplit(args[i],"="))
  if (res[1] == "meta") meta = res[2]
  if (res[1] == "eqtl_normial") eqtl_normial = res[2]
  if (res[1] == "eqtl_permute") eqtl_permute = res[2]
  #if (res[1] == "gwas_sumstat") gwas_sumstat = res[2]   
}

# ## setwd("/sigma4/projects/rnaQTL_hg38/work/85/ce5b0bb1e7cbd639567afda3460949")
# eqtl_normial = "nominal_pass_all_chrom.tsv.gz"
# eqtl_permute = "recount.tsv"
# gwas_sumstat = "T1D.tsv.gz"
# meta = "meta_sumstat.csv"



do_coloc <- function(gwas_sumstat, case_prop, qtl_all, qtl_hit){
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
        input = as.data.frame(input)
        if(nrow(input) > 0){
        dataset1=list(pvalues=input$p_value, type="cc", s = case_prop, N=nrow(gwas))
        dataset2=list(pvalues=input$eqtl_pvalue, type="quant", N=nrow(eqtl))
        result <- coloc.abf( dataset1, dataset2, MAF=input$effect_allele_frequency)
        #qtl_hit$PP.H3.abf[i] = result$summary[5]
        input$ID = input$V1
        input$rsID = input$rs_id
        input$gwas_pvalue = gwas$p_value
        pick = c("ID", "rsID", "eqtl_pvalue", "gwas_pvalue" )
        input_res = input[,pick]
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

## QTLtools
# qtl_all = fread(eqtl_normial, sep = " ")
# qtl_all$rs_id = qtl_all$V8 
# qtl_all$eqtl_pvalue =  qtl_all$V12

# fastQTL
qtl_all = fread(eqtl_normial)
qtl_all$rs_id = qtl_all$V2
qtl_all$eqtl_pvalue =  qtl_all$V7

qtl_all = qtl_all[!is.na(qtl_all$eqtl_pvalue),]
qtl_hit = fread(eqtl_permute)
meta_info = fread(meta)

for( j in 1:nrow(meta_info)){
    gwas_sumstat = as.character(meta_info[j,1])
    case_prop = as.numeric(meta_info[j,2])
    if(!file.exists(gwas_sumstat)){
        cat( paste0(gwas_sumstat, " does not exist! Please check again!\n"))
    }else{
        do_coloc(gwas_sumstat, case_prop, qtl_all, qtl_hit)
    }    
}







