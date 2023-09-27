#!/usr/bin/env nextflow
/*
========================================================================================
                cscQTL : consensus-based circRNA QTL analysis - bed input
========================================================================================
 Circular RNA QTL Analysis Pipeline.
 https://github.com/datngu/cscQTL
 Author: Dat T Nguyen
 Contact: ndat<at>utexas.edu
----------------------------------------------------------------------------------------
*/


/*
To provide highest level of fexlibility for cscQTL user. We provide this a secondary implemtation to run circRNA QTL analysis with given circRNA expression level input in bed formated file.
This assumed that you input 1 bed file for 1 sample. The bed file is 1-based coordinate, with circRNA expression is normalized in CPM and the file name should be the sample_id.bed, otherwise, the program will result an error.
*/

/*
 Define the default parameters
*/ 
params.genotype        = "$baseDir/data/genotype.vcf.gz"
params.meta            = "$baseDir/data/meta.csv"
params.sumstat_files   = "$baseDir/data/sumstat/*.gz"
params.meta_sumstat    = "$baseDir/data/meta_sumstat.csv"
// bed input
params.bed_files       = "$baseDir/data/bed/*.bed"
// out dirs
params.trace_dir       = "trace_dir"
params.outdir          = "results"

// running options
params.chrom           = 1..22 
params.peer            = 1..20 
params.genotype_PCs    = 4 
params.bsj_filter      = 2
params.exp_prop        = 0.3
params.fdr             = 0.05
params.fastqtl_window  = 1000000

// pipeline options
params.coloc           = false





log.info """\
================================================================
                nf-rnaQTL - bed input 
================================================================
    genotype            : $params.genotype 
    meta                : $params.meta
    meta_sumstat        : $params.meta_sumstat
    sumstat_files       : $params.sumstat_files
    bed_files           : $params.bed_files
    outdir              : $params.outdir
    chrom               : $params.chrom
    bsj_filter          : $params.bsj_filter
    exp_prop            : $params.exp_prop
    peer                : $params.peer
    fdr                 : $params.fdr
    fastqtl_window      : $params.fastqtl_window
    genotype_PCs        : $params.genotype_PCs
    coloc               : $params.coloc
================================================================
"""


nextflow.enable.dsl=2

include { VCF_filtering; PCA_genotype } from './main'

workflow {
    // general processing
    chrom_list_ch = channel.from(params.chrom)
    peer_list_ch = channel.from(params.peer)
    bed_ch = channel.fromPath( params.bed_files, checkIfExists: true )

    VCF_filtering(params.genotype, params.meta)
    PCA_genotype(VCF_filtering.out, params.genotype_PCs)

    // data processing and quantile normialization
    BED_merge_samples(params.bsj_filter, params.exp_prop, bed_ch.collect())
    BED_quantile_norm(BED_merge_samples.out, params.meta)
    // preprocessing and circQTL mapping
    BED_covariate_processing(params.meta, BED_quantile_norm.out, PCA_genotype.out, peer_list_ch)
    BED_chrom_splitting(chrom_list_ch, VCF_filtering.out, BED_quantile_norm.out)
    BED_qtl_map_input_ch = BED_covariate_processing.out.combine(BED_chrom_splitting.out)
    BED_qtl_mapping(BED_qtl_map_input_ch, params.fastqtl_window)
    BED_merge_qtl_results(BED_qtl_mapping.out.groupTuple())
    BED_apply_qvalue(params.fdr, BED_merge_qtl_results.out.collect())
    BED_export_all_peer_covariates(BED_covariate_processing.out)
    BED_export_optimal_peer_covariates(BED_apply_qvalue.out, BED_export_all_peer_covariates.out.collect())
    BED_qtl_mapping_nominal(BED_export_optimal_peer_covariates.out, BED_chrom_splitting.out, params.fastqtl_window)
    BED_merge_qtl_results_nominal(BED_qtl_mapping_nominal.out.collect())


    if( params.coloc ){
        sumstat_ch = channel.fromPath( params.sumstat_files, checkIfExists: true )
        COLOC_BED(params.meta_sumstat, BED_apply_qvalue.out, BED_merge_qtl_results_nominal.out, sumstat_ch.collect())
        
    }
}


// steps
process COLOC_BED {
    container 'ndatth/qtl-package:v0.0.0'
    publishDir "${params.outdir}/coloc_cscQTL_bed", mode: 'copy', overwrite: true
    memory '32 GB'
    cpus 8

    input:
    path meta_sumstat
    path eqtl_normial
    path eqtl_permute_fdr
    path sumstat_files

    output:
    path "coloc*"

    script:
    """
    coloc.R meta=$meta_sumstat eqtl_normial=nominal_pass_all_chrom.tsv.gz eqtl_permute=cscQTL_bed.tsv
    """
}


process BED_merge_samples {
    container 'ndatth/rna-tools:v0.0.0'
    publishDir "${params.outdir}/cscQTL_bed", mode: 'copy', overwrite: true
    memory '8 GB'

    input:
    val bsj_filter
    val exp_prop
    path beds

    output:
    path "bed_merged_bsj_${bsj_filter}.tsv"

    script:
    """
    bed_merge_results.R "bed_merged_bsj_${bsj_filter}.tsv" ${bsj_filter} ${exp_prop} ${beds}
    """
}


process BED_quantile_norm { 
    container 'ndatth/rna-tools:v0.0.0'
    publishDir "${params.outdir}/cscQTL_bed", mode: 'copy', overwrite: true
    memory '8 GB'
    
    input:
    path "res_merged.tsv"
    path "meta.csv"
    
    output:
    path "bed_quantile_norm.tsv"
 
    script:
    """
    quantile_norm.R res_merged.tsv meta.csv bed_quantile_norm.tsv
    """
}



process BED_covariate_processing { 
    
    publishDir "${params.outdir}/cscQTL_bed_qtl_input", mode: 'copy', overwrite: true
    container 'ndatth/qtl-package:v0.0.0'
    memory '8 GB'

    input:
    path meta
    path exp_quatile_normalized
    path genotype_PCs 
    val n
 
    output:
    tuple val("${n}"), path("${n}_peer_covariates.tsv")
 
    script:
    """
    covariate_processing.R ${meta} ${exp_quatile_normalized} ${genotype_PCs} ${n} ${n}_peer_covariates.tsv
    """
}




process BED_chrom_splitting {

    publishDir "${params.outdir}/cscQTL_bed_qtl_input", mode: 'copy', overwrite: true
    container 'ndatth/qtl-package:v0.0.0'
    memory '2 GB'

    input:
    val chr
    tuple path("in.vcf.gz"), path("in.vcf.gz.tbi")
    path exp

    output:
    tuple val("${chr}"), path("${chr}.vcf.gz"), path("${chr}.vcf.gz.tbi"), path("${chr}.bed.gz"), path("${chr}.bed.gz.tbi")
 
    script:
    """
    # vcf files
    bcftools view in.vcf.gz --regions ${chr} -Oz -o ${chr}.vcf.gz
    bcftools index -t ${chr}.vcf.gz
    # exp files
    grep -w "^#Chr" ${exp} > ${chr}.bed
    grep -w ^${chr} ${exp} | sort -nk2 >> ${chr}.bed
    cat ${chr}.bed | bgzip > ${chr}.bed.gz
    tabix -f ${chr}.bed.gz
    rm ${chr}.bed
    """
}



process BED_qtl_mapping {

    publishDir "${params.trace_dir}/cscQTL_bed_qtl_mapping", mode: 'symlink', overwrite: true
    container 'ndatth/qtl-package:v0.0.0'
    memory '2 GB'

    input:
    tuple val("n"), path("covariates.tsv"), val("chr"), path("chr.vcf.gz"), path("chr.vcf.gz.tbi"), path("chr.bed.gz"), path("chr.bed.gz.tbi")
    val window_size

    output:
    
    tuple val("${n}"), path("peer_${n}_chr_${chr}_permute_pass.tsv.gz")
 
    script:
    """
    fastQTL --vcf chr.vcf.gz --bed chr.bed.gz --cov covariates.tsv --window ${window_size} --permute 1000 10000 --chunk 1 1 --out "peer_${n}_chr_${chr}_permute_pass.tsv"

    gzip -f peer_${n}_chr_${chr}_permute_pass.tsv
    """
}



process BED_merge_qtl_results {

    publishDir "${params.trace_dir}/cscQTL_bed_qtl_mapping_merged", mode: 'symlink', overwrite: true
    container 'ndatth/qtl-package:v0.0.0'
    memory '8 GB'

    input:
    tuple val("n"), path(qtl_res)

    output:
    
    path ("${n}_peer_permute_pass_all_chrom.tsv.gz")
 
    script:
    """
    zcat ${qtl_res} > ${n}_peer_permute_pass_all_chrom.tsv

    gzip -f ${n}_peer_permute_pass_all_chrom.tsv
    """
}


process BED_apply_qvalue {

    publishDir "${params.outdir}/qtl_mapping_apply_qvalue", mode: 'copy', overwrite: true
    container 'ndatth/qtl-package:v0.0.0'
    memory '8 GB'

    input:
    val fdr
    path merge_qtl_results

    output:
    
    tuple path("cscQTL_bed"), path("cscQTL_bed.tsv"), path("cscQTL_bed.stat"), path("cscQTL_bed.png"), path("cscQTL_bed_raw.tsv") 
 
    script:
    """
    fastQTL_apply_qvalue.R cscQTL_bed ${fdr} ${merge_qtl_results}
    """
}


process BED_export_all_peer_covariates { 
    
    publishDir "${params.trace_dir}/cscQTL_bed_qtl_covariates", mode: 'symlink', overwrite: true
    container 'ndatth/qtl-package:v0.0.0'
    memory '4 GB'

    input:
    tuple val(n), path(n_peer_covariates)
 
    output:
    path "${n}_peer_covariates_exported.tsv"
 
    script:
    """
    cat ${n_peer_covariates} > ${n}_peer_covariates_exported.tsv
    """
}



process BED_export_optimal_peer_covariates {

    publishDir "${params.trace_dir}/cscQTL_bed_qtl_covariates", mode: 'symlink', overwrite: true
    container 'ndatth/qtl-package:v0.0.0'
    memory '2 GB'

    input:
    tuple path("cscQTL_bed"), path("cscQTL_bed.tsv"), path("cscQTL_bed.stat"), path("cscQTL_bed.png"), path("cscQTL_bed_raw.tsv") 
    path all_peer_factor_covariates

    output:
    
    path("optimal_covariates.tsv") 
 
    script:
    """
    n="\$(cat cscQTL_bed)"
    cat \${n}_peer_covariates_exported.tsv > optimal_covariates.tsv
    """
}



process BED_qtl_mapping_nominal {

    publishDir "${params.outdir}/cscQTL_bed_qtl_mapping_nominal", mode: 'copy', overwrite: true
    container 'ndatth/qtl-package:v0.0.0'
    memory '2 GB'

    input:
    path("covariates.tsv")
    tuple val("chr"), path("chr.vcf.gz"), path("chr.vcf.gz.tbi"), path("chr.bed.gz"), path("chr.bed.gz.tbi")
    val window_size

    output:
    
    path("chr_${chr}_nominal_pass_best_peer_factors.tsv.gz")
 
    script:
    """
    zcat chr.bed.gz | awk '{ \$4=\$4" . +"; print \$0 }' | tr " " "\t" | bgzip -c > QTLtools_chr.bed.gz
    tabix -f QTLtools_chr.bed.gz
    QTLtools cis --vcf chr.vcf.gz --bed QTLtools_chr.bed.gz --cov covariates.tsv --window ${window_size} --nominal 0.1 --chunk 1 1 --out chr_${chr}_nominal_pass_best_peer_factors.tsv
    gzip -f chr_${chr}_nominal_pass_best_peer_factors.tsv
    """
}




process BED_merge_qtl_results_nominal {

    publishDir "${params.trace_dir}/cscQTL_bed_qtl_mapping_nominal", mode: 'symlink', overwrite: true
    container 'ndatth/qtl-package:v0.0.0'
    memory '8 GB'

    input:
    path(qtl_nominal_res)

    output:
    
    path ("nominal_pass_all_chrom.tsv.gz")
 
    script:
    """
    zcat ${qtl_nominal_res} > nominal_pass_all_chrom.tsv
    gzip -f nominal_pass_all_chrom.tsv
    """
}

