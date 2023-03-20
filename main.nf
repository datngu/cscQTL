#!/usr/bin/env nextflow
/*
========================================================================================
                    cscQTL : consensus-based circRNA QTL analysis
========================================================================================
 Circular RNA QTL Analysis Pipeline.
 https://github.com/datngu/cscQTL
 Author: Dat T Nguyen
 Contact: ndat<at>utexas.edu
----------------------------------------------------------------------------------------
*/


/*
 Define the default parameters
*/ 
params.genome          = "$baseDir/data/ref/genome.fa"
params.cdna            = "$baseDir/data/ref/cdna.fa"
params.bsj             = "$baseDir/data/ref/bsj.fa"
params.annotation      = "$baseDir/data/ref/annotation.gtf"
params.reads           = "$baseDir/data/reads/*_{1,2}.fastq.gz"
params.genotype        = "$baseDir/data/genotype.vcf.gz"
params.meta            = "$baseDir/data/meta.csv"
params.sumstat_files   = "$baseDir/data/sumstat/*.gz"
params.meta_sumstat    = "$baseDir/data/meta_sumstat.csv"
params.trace_dir       = "trace_dir"
params.outdir          = "results"

// running options
params.chrom           = 1..22 
params.peer            = 1..20 
params.genotype_PCs    = 4 
params.bsj_filter      = 2
params.consensus       = 1
params.exp_prop        = 0.3
params.maf             = 0.05
params.fdr             = 0.05
params.fastqtl_window  = 1000000

// pipeline options
params.cscqtl          = true
params.coloc           = false
params.circall         = false
params.ciri2           = false
params.circexplorer2   = false
params.salmon          = false



log.info """\
================================================================
                      nf-rnaQTL 
================================================================
    genome              : $params.genome
    cdna                : $params.cdna
    bsj                 : $params.bsj
    annotation          : $params.annotation
    reads               : $params.reads
    genotype            : $params.genotype 
    meta                : $params.meta
    meta_sumstat        : $params.meta_sumstat
    sumstat_files       : $params.sumstat_files
    outdir              : $params.outdir
    chrom               : $params.chrom
    bsj_filter          : $params.bsj_filter
    consensus           : $params.consensus
    exp_prop            : $params.exp_prop
    peer                : $params.peer
    maf                 : $params.maf
    fdr                 : $params.fdr
    fastqtl_window      : $params.fastqtl_window
    genotype_PCs        : $params.genotype_PCs
    circall             : $params.circall
    ciri2               : $params.ciri2
    circexplorer2       : $params.circexplorer2
    cscqtl              : $params.cscqtl
    salmon              : $params.salmon
    coloc               : $params.coloc
================================================================
"""


nextflow.enable.dsl=2

include { RECOUNT_get_consensus; RENAME_circall; RENAME_ciri2; RENAME_circexp2; RECOUNT_get_peseudo_seqs; RECOUNT_index_peseudo_seqs; RECOUNT_mapping; RECOUNT_count; RECOUNT_merge_samples; RECOUNT_quantile_norm; RECOUNT_covariate_processing; RECOUNT_chrom_splitting; RECOUNT_qtl_mapping; RECOUNT_merge_qtl_results; RECOUNT_apply_qvalue; RECOUNT_export_all_peer_covariates; RECOUNT_export_optimal_peer_covariates; RECOUNT_qtl_mapping_nominal; RECOUNT_merge_qtl_results_nominal } from './module/recount'


workflow {
    // general processing
    read_pairs_ch = channel.fromFilePairs( params.reads, checkIfExists: true )
    chrom_list_ch = channel.from(params.chrom)
    peer_list_ch = channel.from(params.peer)
    VCF_filtering(params.genotype, params.meta)
    PCA_genotype(VCF_filtering.out, params.genotype_PCs)

    // run cscQTL - the functional implemented is recount - please noted for debugging
    if( params.cscqtl ){
        // running Circall pipeline
        CIRCALL_index_cdna(params.cdna)
        CIRCALL_index_bsj(params.bsj)
        CIRCALL_sqlite_ano(params.annotation)
        CIRCALL_pipeline(params.genome, params.cdna, CIRCALL_index_cdna.out, CIRCALL_index_bsj.out, CIRCALL_sqlite_ano.out, read_pairs_ch)
        // running CIRI2 pipeline
        BWA_index_genome(params.genome)
        CIRI2_pipeline(params.genome, params.annotation, BWA_index_genome.out, read_pairs_ch)
        // running CE2
        STAR_index_genome(params.genome, params.annotation)
        STAR_mapping(STAR_index_genome.out, read_pairs_ch)   
        CIRCEXP2_generate_annotation(params.annotation)
        CIRCEXP2_pipeline(params.genome, CIRCEXP2_generate_annotation.out, STAR_mapping.out)

        // rename and merging input of the 3 pipelines
        RENAME_circall(CIRCALL_pipeline.out)
        RENAME_ciri2(CIRI2_pipeline.out)
        RENAME_circexp2(CIRCEXP2_pipeline.out)
        RECOUNT_get_consensus(params.annotation, RENAME_circall.out.collect(), RENAME_ciri2.out.collect(), RENAME_circexp2.out.collect())
        // generating peseudo references
        RECOUNT_get_peseudo_seqs(params.genome, RECOUNT_get_consensus.out)
        RECOUNT_index_peseudo_seqs(RECOUNT_get_peseudo_seqs.out)
        // mapping and counting BSJ reads
        RECOUNT_mapping(CIRCALL_index_cdna.out, RECOUNT_index_peseudo_seqs.out, read_pairs_ch)
        RECOUNT_count(RECOUNT_mapping.out)
        // data processing and quantile normialization
        RECOUNT_merge_samples(params.bsj_filter, params.exp_prop, RECOUNT_count.out.collect())
        RECOUNT_quantile_norm(RECOUNT_merge_samples.out, params.meta)
        // preprocessing and circQTL mapping
        RECOUNT_covariate_processing(params.meta, RECOUNT_quantile_norm.out, PCA_genotype.out, peer_list_ch)
        RECOUNT_chrom_splitting(chrom_list_ch, VCF_filtering.out, RECOUNT_quantile_norm.out)
        RECOUNT_qtl_map_input_ch = RECOUNT_covariate_processing.out.combine(RECOUNT_chrom_splitting.out)
        RECOUNT_qtl_mapping(RECOUNT_qtl_map_input_ch, params.fastqtl_window)
        RECOUNT_merge_qtl_results(RECOUNT_qtl_mapping.out.groupTuple())
        RECOUNT_apply_qvalue(params.fdr, RECOUNT_merge_qtl_results.out.collect())
        RECOUNT_export_all_peer_covariates(RECOUNT_covariate_processing.out)
        RECOUNT_export_optimal_peer_covariates(RECOUNT_apply_qvalue.out, RECOUNT_export_all_peer_covariates.out.collect())
        RECOUNT_qtl_mapping_nominal(RECOUNT_export_optimal_peer_covariates.out, RECOUNT_chrom_splitting.out, params.fastqtl_window)
        RECOUNT_merge_qtl_results_nominal(RECOUNT_qtl_mapping_nominal.out.collect())
    }



    // This step is shared bw single methods!
    if( params.circall || params.ciri2 || params.circexplorer2){
        LIBSIZE_count(read_pairs_ch)
        LIBSIZE_merge_samples(LIBSIZE_count.out.collect())
    }


    // running circQLT for Circall
    if( params.circall ){
        CIRCALL_merge_samples(LIBSIZE_merge_samples.out, params.bsj_filter, params.exp_prop, CIRCALL_pipeline.out.collect())
        CIRCALL_quantile_norm(CIRCALL_merge_samples.out, params.meta)
        //
        CIRCALL_covariate_processing(params.meta, CIRCALL_quantile_norm.out, PCA_genotype.out, peer_list_ch)
        CIRCALL_chrom_splitting(chrom_list_ch, VCF_filtering.out, CIRCALL_quantile_norm.out)
        CIRCALL_qtl_map_input_ch = CIRCALL_covariate_processing.out.combine(CIRCALL_chrom_splitting.out)
        CIRCALL_qtl_mapping(CIRCALL_qtl_map_input_ch, params.fastqtl_window)
        CIRCALL_merge_qtl_results(CIRCALL_qtl_mapping.out.groupTuple())
        //
        CIRCALL_apply_qvalue(params.fdr, CIRCALL_merge_qtl_results.out.collect())
        CIRCALL_export_all_peer_covariates(CIRCALL_covariate_processing.out)
        CIRCALL_export_optimal_peer_covariates(CIRCALL_apply_qvalue.out, CIRCALL_export_all_peer_covariates.out.collect())
        CIRCALL_qtl_mapping_nominal(CIRCALL_export_optimal_peer_covariates.out, CIRCALL_chrom_splitting.out, params.fastqtl_window)
        CIRCALL_merge_qtl_results_nominal(CIRCALL_qtl_mapping_nominal.out.collect())
    }

    
    // running circQLT for CIRI2
    if( params.ciri2 ){
        CIRI2_merge_samples(LIBSIZE_merge_samples.out, params.bsj_filter, params.exp_prop, CIRI2_pipeline.out.collect())
        CIRI2_quantile_norm(CIRI2_merge_samples.out, params.meta)
        //
        CIRI2_covariate_processing(params.meta, CIRI2_quantile_norm.out, PCA_genotype.out, peer_list_ch)
        CIRI2_chrom_splitting(chrom_list_ch, VCF_filtering.out, CIRI2_quantile_norm.out)
        CIRI2_qtl_map_input_ch = CIRI2_covariate_processing.out.combine(CIRI2_chrom_splitting.out)
        CIRI2_qtl_mapping(CIRI2_qtl_map_input_ch, params.fastqtl_window)
        CIRI2_merge_qtl_results(CIRI2_qtl_mapping.out.groupTuple())
        //
        CIRI2_apply_qvalue(params.fdr, CIRI2_merge_qtl_results.out.collect())
        CIRI2_export_all_peer_covariates(CIRI2_covariate_processing.out)
        CIRI2_export_optimal_peer_covariates(CIRI2_apply_qvalue.out, CIRI2_export_all_peer_covariates.out.collect())
        CIRI2_qtl_mapping_nominal(CIRI2_export_optimal_peer_covariates.out, CIRI2_chrom_splitting.out, params.fastqtl_window)
        CIRI2_merge_qtl_results_nominal(CIRI2_qtl_mapping_nominal.out.collect())               
    }

    

    // running circQLT for circexplorer2
    if( params.circexplorer2 ){
        CIRCEXP2_merge_samples(LIBSIZE_merge_samples.out, params.bsj_filter, params.exp_prop, CIRCEXP2_pipeline.out.collect())
        CIRCEXP2_quantile_norm(CIRCEXP2_merge_samples.out, params.meta)
        //
        CIRCEXP2_covariate_processing(params.meta, CIRCEXP2_quantile_norm.out, PCA_genotype.out, peer_list_ch)
        CIRCEXP2_chrom_splitting(chrom_list_ch, VCF_filtering.out, CIRCEXP2_quantile_norm.out)
        CIRCEXP2_qtl_map_input_ch = CIRCEXP2_covariate_processing.out.combine(CIRCEXP2_chrom_splitting.out)
        CIRCEXP2_qtl_mapping(CIRCEXP2_qtl_map_input_ch, params.fastqtl_window)
        CIRCEXP2_merge_qtl_results(CIRCEXP2_qtl_mapping.out.groupTuple())
        //
        CIRCEXP2_apply_qvalue(params.fdr, CIRCEXP2_merge_qtl_results.out.collect())
        CIRCEXP2_export_all_peer_covariates(CIRCEXP2_covariate_processing.out)
        CIRCEXP2_export_optimal_peer_covariates(CIRCEXP2_apply_qvalue.out, CIRCEXP2_export_all_peer_covariates.out.collect())
        CIRCEXP2_qtl_mapping_nominal(CIRCEXP2_export_optimal_peer_covariates.out, CIRCEXP2_chrom_splitting.out, params.fastqtl_window)
        CIRCEXP2_merge_qtl_results_nominal(CIRCEXP2_qtl_mapping_nominal.out.collect())
    }



    if( params.salmon ){
        SALMON_index_cdna(params.cdna)
        SALMON_quantify(SALMON_index_cdna.out, read_pairs_ch)
        SALMON_merge_samples(params.annotation, SALMON_quantify.out.collect())
        SALMON_quantile_norm(SALMON_merge_samples.out, params.meta, params.exp_prop, params.annotation)
        //
        SALMON_covariate_processing(params.meta, SALMON_quantile_norm.out, PCA_genotype.out, peer_list_ch)
        SALMON_chrom_splitting(chrom_list_ch, VCF_filtering.out, SALMON_quantile_norm.out)
        SALMON_qtl_map_input_ch = SALMON_covariate_processing.out.combine(SALMON_chrom_splitting.out)
        SALMON_qtl_mapping(SALMON_qtl_map_input_ch, params.fastqtl_window)
        SALMON_merge_qtl_results(SALMON_qtl_mapping.out.groupTuple())
        SALMON_apply_qvalue(params.fdr, SALMON_merge_qtl_results.out.collect())
        SALMON_export_all_peer_covariates(SALMON_covariate_processing.out)
        SALMON_export_optimal_peer_covariates(SALMON_apply_qvalue.out, SALMON_export_all_peer_covariates.out.collect())
        SALMON_qtl_mapping_nominal(SALMON_export_optimal_peer_covariates.out, SALMON_chrom_splitting.out, params.fastqtl_window)
        SALMON_merge_qtl_results_nominal(SALMON_qtl_mapping_nominal.out.collect())
    }

    if( params.coloc ){
        sumstat_ch = channel.fromPath( params.sumstat_files, checkIfExists: true )
        
        if( params.cscqtl ){
            COLOC_recount(params.meta_sumstat, RECOUNT_apply_qvalue.out, RECOUNT_merge_qtl_results_nominal.out, sumstat_ch.collect()) 
        }
        
    }

}

//////////////////////// done pipeline!


process COLOC_circall {
    container 'ndatth/qtl-package:v0.0.0'
    publishDir "${params.outdir}/coloc_circall", mode: 'copy', overwrite: true
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
    coloc.R meta=$meta_sumstat eqtl_normial=nominal_pass_all_chrom.tsv.gz eqtl_permute=circall.tsv
    """
}

process COLOC_ciri2 {
    container 'ndatth/qtl-package:v0.0.0'
    publishDir "${params.outdir}/coloc_ciri2", mode: 'copy', overwrite: true
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
    coloc.R meta=$meta_sumstat eqtl_normial=nominal_pass_all_chrom.tsv.gz eqtl_permute=ciri2.tsv
    """
}

process COLOC_circexp2 {
    container 'ndatth/qtl-package:v0.0.0'
    publishDir "${params.outdir}/coloc_circexp2", mode: 'copy', overwrite: true
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
    coloc.R meta=$meta_sumstat eqtl_normial=nominal_pass_all_chrom.tsv.gz eqtl_permute=circexp2.tsv
    """
}


process COLOC_recount {
    container 'ndatth/qtl-package:v0.0.0'
    publishDir "${params.outdir}/coloc_recount", mode: 'copy', overwrite: true
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
    coloc.R meta=$meta_sumstat eqtl_normial=nominal_pass_all_chrom.tsv.gz eqtl_permute=recount.tsv
    """
}



process COLOC_salmon {

    container 'ndatth/qtl-package:v0.0.0'
    publishDir "${params.outdir}/coloc_salmon", mode: 'copy', overwrite: true
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
    coloc.R meta=$meta_sumstat eqtl_normial=nominal_pass_all_chrom.tsv.gz eqtl_permute=salmon.tsv
    """
}





process LIBSIZE_count {

    container 'ndatth/rna-tools:v0.0.0'
    publishDir "${params.trace_dir}/libsize", mode: 'symlink', overwrite: true

    input:
    tuple val(pair_id), path(reads)
 
    output:
    path "${pair_id}_libsize.txt"
 
    script:
    """
    libsize_count.py ${pair_id} ${reads[0]} ${pair_id}_libsize.txt
    """
}


process LIBSIZE_merge_samples {

    container 'ndatth/rna-tools:v0.0.0'
    publishDir "${params.trace_dir}/libsize", mode: 'symlink', overwrite: true

    input:
    path libsize_count_outputs
 
    output:
    path "libsize.tsv"
 
    script:
    """
    cat ${libsize_count_outputs} > libsize.tsv
    """
}



// GENOTYPE PROCESSING
process VCF_filtering { 
    
    publishDir "${params.trace_dir}/vcf_filtering", mode: 'symlink', overwrite: true
    container 'ndatth/qtl-package:v0.0.0'
    memory '8 GB'
    
    input:
    path "in.vcf"
    path "meta.csv"
 
    output:
    tuple path("genotype_filtered.vcf.gz"), path("genotype_filtered.vcf.gz.tbi")
 
    script:
    """
    #"in.vcf"=tcel.vcf.gz
    #meta=/sigma4/projects/nf-circall-qtl/data/meta.csv
    plink --vcf in.vcf --make-bed --out genotype_raw --memory 8000
    plink --bfile genotype_raw --maf 0.05 --hwe 1e-6 --memory 8000 --make-bed --const-fid --out genotype_QCed
    plink --bfile genotype_QCed --memory 8000 --recode vcf-fid --out genotype_filtered_tem
    grep -v ^rna_id meta.csv | cut -d , -f 2 > genotype_sample.txt
    ## filtering sample and remove chr (if needed)
    bcftools view genotype_filtered_tem.vcf -S genotype_sample.txt | sed 's/chr//g' | bgzip > genotype_filtered.vcf.gz
    bcftools index -t genotype_filtered.vcf.gz
    rm genotype_raw*
    rm genotype_QCed*
    rm genotype_filtered_tem*
    """
}


process PCA_genotype { 
    
    publishDir "${params.trace_dir}/vcf_filtering", mode: 'symlink', overwrite: true
    container 'ndatth/qtl-package:v0.0.0'
    memory '8 GB'
    
    input:
    tuple path("in.vcf.gz"), path("in.vcf.gz.tbi")
    val "genotype_PCs"

    output:
    path "genotype_PCs.tsv"
 
    script:
    """
    mkdir plink
    plink --vcf in.vcf.gz \
      --vcf-half-call 'haploid' \
      --allow-extra-chr \
      --make-bed --const-fid \
      --out plink/tem \
      --threads 1 \
      --memory 8000


    # prunning
    plink \
      --bfile plink/tem \
      --allow-extra-chr \
      --indep-pairwise 200 50 0.25 \
      --memory 8000 \
      --out plink/tem_pruned


    # calculate the first n PCs
    plink \
      --bfile plink/tem \
      --allow-extra-chr \
      --extract plink/tem_pruned.prune.in \
      --memory 8000 \
      --pca ${genotype_PCs} \
      --out plink/genotype_PCs

    cp plink/genotype_PCs.eigenvec genotype_PCs.tsv
    rm -rf plink
    """
}



///////////////////////////////////////////////////////////////
// CIRCALL
process CIRCALL_index_cdna {
    container 'ndatth/rna-tools:v0.0.0'

    input:
    path "cdna.fa"

    output:
    path "index_cdna"

    script:
    """
    TxIndexer -t cdna.fa -o index_cdna
    """
}

process CIRCALL_index_bsj {
    container 'ndatth/rna-tools:v0.0.0'

    input:
    path "bsj.fa"

    output:
    path "index_bsj"

    script:
    """
    TxIndexer -t bsj.fa -o index_bsj
    """
}


process CIRCALL_sqlite_ano {
    container 'ndatth/rna-tools:v0.0.0'

    input:
    path "annotation.gtf"

    output:
    path "annotation.sqlite"

    script:
    """
    createSqlite.R annotation.gtf annotation.sqlite
    """
}



process CIRCALL_pipeline {
    container 'ndatth/rna-tools:v0.0.0'
    publishDir "${params.outdir}/circall" , mode: 'copy', overwrite: true
    cpus 32
    memory '45 GB'

    input:
    path "genome.fa"
    path "cdna.fa"
    path "index_cdna"
    path "index_bsj"
    path "annotation.sqlite"
    tuple val(pair_id), path(reads)
 
    output:
    path "${pair_id}"
 
    script:
    """
    Circall.sh -genome genome.fa -gtfSqlite annotation.sqlite -txFasta cdna.fa -txIdx index_cdna -bsjIdx index_bsj -read1 ${reads[0]} -read2 ${reads[1]} -o ${pair_id}_dir -p ${task.cpus} -td TRUE -dep Circall/Data/Circall_depdata_human.RData
    mv ${pair_id}_dir/Sample_Circall_final.txt ${pair_id}
    """
}


process CIRCALL_merge_samples { 
    container 'ndatth/rna-tools:v0.0.0'
    publishDir "${params.outdir}/circall", mode: 'copy', overwrite: true

    input:
    path "libsize.tsv"
    val bsj_filter
    val exp_prop
    path circall_outputs
 
    output:
    path "circall_res_merged_bsj_${bsj_filter}.tsv"

    script:
    """
    circall_merge_results.R circall_res_merged_bsj_${bsj_filter}.tsv libsize.tsv ${bsj_filter} ${exp_prop} ${circall_outputs}
    """
}



process CIRCALL_quantile_norm { 
    container 'ndatth/rna-tools:v0.0.0'
    publishDir "${params.outdir}/quantile_norm", mode: 'copy', overwrite: true
    memory '8 GB'
    
    input:
    path "res_merged.tsv"
    path "meta.csv"
 
    output:
    path "circall_quantile_norm.tsv"
 
    script:
    """
    quantile_norm.R res_merged.tsv meta.csv circall_quantile_norm.tsv
    """
}


process CIRCALL_covariate_processing { 
    
    publishDir "${params.trace_dir}/circall_qtl_input", mode: 'symlink', overwrite: true
    container 'ndatth/qtl-package:v0.0.0'
    memory '4 GB'

    input:
    path meta
    path exp_quatile_normalized
    path genotype_PCs 
    val n
 
    output:
    tuple val("${n}"), path("${n}_peer_covariates.tsv")
 
    script:
    """
    ## covariate_processing.R meta.csv exp_quatile_normalized.tsv genotype_PCs.tsv number_hidden_factors out_fn
    covariate_processing.R ${meta} ${exp_quatile_normalized} ${genotype_PCs} ${n} ${n}_peer_covariates.tsv
    """
}


process CIRCALL_chrom_splitting {

    publishDir "${params.outdir}/circall_qtl_input", mode: 'copy', overwrite: true
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



process CIRCALL_qtl_mapping {

    publishDir "${params.trace_dir}/circall_qtl_mapping", mode: 'symlink', overwrite: true
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


process CIRCALL_merge_qtl_results {

    publishDir "${params.trace_dir}/circall_qtl_mapping_merged", mode: 'symlink', overwrite: true
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


process CIRCALL_apply_qvalue {

    publishDir "${params.outdir}/qtl_mapping_apply_qvalue", mode: 'copy', overwrite: true
    container 'ndatth/qtl-package:v0.0.0'
    memory '8 GB'

    input:
    val fdr
    path merge_qtl_results

    output:
    
    tuple path("circall"), path("circall.tsv"), path("circall.stat"), path("circall.png"), path("circall_raw.tsv")
 
    script:
    """
    fastQTL_apply_qvalue.R circall ${fdr} ${merge_qtl_results}
    """
}

process CIRCALL_export_all_peer_covariates { 
    
    publishDir "${params.trace_dir}/circall_qtl_covariates", mode: 'symlink', overwrite: true
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


process CIRCALL_export_optimal_peer_covariates {

    publishDir "${params.trace_dir}/circall_qtl_covariates", mode: 'symlink', overwrite: true
    container 'ndatth/qtl-package:v0.0.0'
    memory '2 GB'

    input:
    tuple path("circall"), path("circall.tsv"), path("circall.stat"), path("circall.png"), path("circall_raw.tsv")
    path all_peer_factor_covariates

    output:
    
    path("optimal_covariates.tsv") 
 
    script:
    """
    n="\$(cat circall)"
    cat \${n}_peer_covariates_exported.tsv > optimal_covariates.tsv
    """
}


process CIRCALL_qtl_mapping_nominal {

    publishDir "${params.outdir}/circall_qtl_mapping_nominal", mode: 'copy', overwrite: true
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


process CIRCALL_merge_qtl_results_nominal {

    publishDir "${params.trace_dir}/circall_qtl_mapping_nominal", mode: 'symlink', overwrite: true
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


///////////////////////////////////////////////////////////////
// CIRI2


process BWA_index_genome {
    container 'ndatth/rna-tools:v0.0.0'
    publishDir "${params.trace_dir}/bwa_index", mode: 'symlink', overwrite: true
    memory '8 GB'

    input:
    path "genome.fa"

    output:
    path "genome.fa*"

    script:
    """
    bwa index -a bwtsw genome.fa
    """
}


process CIRI2_pipeline {
    container 'ndatth/rna-tools:v0.0.0'
    publishDir "${params.outdir}/ciri2", mode: 'copy', overwrite: true
    cpus 32
    memory '45 GB'

    input:
    path "genome.fa"
    path "annotation.gtf"
    path bwa_index
    tuple val(pair_id), path(reads)

    output:
    path "${pair_id}"

    script:
    """
    bwa mem -t ${task.cpus} -T 19 genome.fa ${reads[0]} ${reads[1]} 1> ${pair_id}.sam 2> ${pair_id}.log
    ## cpu number of ciri2 divide by 2 for RAM requirements
    perl /bin/CIRI2.pl -I ${pair_id}.sam -O ${pair_id} -F genome.fa -A annotation.gtf -T 4
    rm ${pair_id}.sam
    """
}


process CIRI2_merge_samples { 
    container 'ndatth/rna-tools:v0.0.0'
    publishDir "${params.outdir}/ciri2", mode: 'copy', overwrite: true

    input:
    path "libsize.tsv"
    val bsj_filter
    val exp_prop
    path ciri2_outputs
 
    output:
    path "ciri2_res_merged_bsj_${bsj_filter}.tsv"

    script:
    """
    ciri2_merge_results.R ciri2_res_merged_bsj_${bsj_filter}.tsv libsize.tsv ${bsj_filter} ${exp_prop} ${ciri2_outputs}
    """
}



process CIRI2_quantile_norm { 
    container 'ndatth/rna-tools:v0.0.0'
    publishDir "${params.outdir}/quantile_norm", mode: 'copy', overwrite: true
    memory '8 GB'
    
    input:
    path "res_merged.tsv"
    path "meta.csv"
    
    output:
    path "ciri2_quantile_norm.tsv"
 
    script:
    """
    quantile_norm.R res_merged.tsv meta.csv ciri2_quantile_norm.tsv
    """
}




process CIRI2_covariate_processing { 
    
    publishDir "${params.trace_dir}/ciri2_qtl_input"
    container 'ndatth/qtl-package:v0.0.0'
    memory '4 GB'

    input:
    path meta
    path exp_quatile_normalized
    path genotype_PCs 
    val n
 
    output:
    tuple val("${n}"), path("${n}_peer_covariates.tsv")
 
    script:
    """
    ## covariate_processing.R meta.csv exp_quatile_normalized.tsv genotype_PCs.tsv number_hidden_factors out_fn
    covariate_processing.R ${meta} ${exp_quatile_normalized} ${genotype_PCs} ${n} ${n}_peer_covariates.tsv
    """
}


process CIRI2_chrom_splitting {

    publishDir "${params.outdir}/ciri2_qtl_input", mode: 'copy', overwrite: true
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



process CIRI2_qtl_mapping {

    publishDir "${params.trace_dir}/ciri2_qtl_mapping", mode: 'symlink', overwrite: true
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


process CIRI2_merge_qtl_results {

    publishDir "${params.trace_dir}/ciri2_qtl_mapping_merged", mode: 'symlink', overwrite: true
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




process CIRI2_apply_qvalue {

    publishDir "${params.outdir}/qtl_mapping_apply_qvalue", mode: 'copy', overwrite: true
    container 'ndatth/qtl-package:v0.0.0'
    memory '8 GB'

    input:
    val fdr
    path merge_qtl_results

    output:
    
    tuple path("ciri2"), path("ciri2.tsv"), path("ciri2.stat"), path("ciri2.png"), path("ciri2_raw.tsv") 
 
    script:
    """
    fastQTL_apply_qvalue.R ciri2 ${fdr} ${merge_qtl_results}
    """
}


process CIRI2_export_all_peer_covariates { 
    
    publishDir "${params.trace_dir}/ciri2_qtl_covariates", mode: 'symlink', overwrite: true
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


process CIRI2_export_optimal_peer_covariates {

    publishDir "${params.trace_dir}/ciri2_qtl_covariates", mode: 'symlink', overwrite: true
    container 'ndatth/qtl-package:v0.0.0'
    memory '2 GB'

    input:
    tuple path("ciri2"), path("ciri2.tsv"), path("ciri2.stat"), path("ciri2.png"), path("ciri2_raw.tsv") 
    path all_peer_factor_covariates

    output:
    
    path("optimal_covariates.tsv") 
 
    script:
    """
    n="\$(cat ciri2)"
    cat \${n}_peer_covariates_exported.tsv > optimal_covariates.tsv
    """
}


process CIRI2_qtl_mapping_nominal {

    publishDir "${params.outdir}/ciri2_qtl_mapping_nominal", mode: 'copy', overwrite: true
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


process CIRI2_merge_qtl_results_nominal {

    publishDir "${params.trace_dir}/ciri2_qtl_mapping_nominal", mode: 'symlink', overwrite: true
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






///////////////////////////////////////////////////////////////
// STAR - (for CIRCexplorer2)

// STAR_index_genome(params.genome, annotation.gtf)
process STAR_index_genome {
    container 'ndatth/rna-tools:v0.0.0'
    publishDir "${params.trace_dir}/star_index", mode: 'symlink', overwrite: true
    cpus 32
    memory '45 GB'
    

    input:
    path "genome.fa"
    path "annotation.gtf"

    output:
    path "STAR_genome_index"

    script:
    """
    mkdir STAR_genome_index
    STAR --runMode genomeGenerate --genomeDir STAR_genome_index --genomeFastaFiles genome.fa --sjdbGTFfile annotation.gtf --sjdbOverhang 100 --runThreadN ${task.cpus}
    """
}


// STAR_mapping(STAR_index_genome.out, read_pairs_ch)   

process STAR_mapping {
    container 'ndatth/rna-tools:v0.0.0'
    publishDir "${params.trace_dir}/star_mapping", mode: 'symlink', overwrite: true
    cpus 32
    memory '45 GB'

    input:
    path "STAR_genome_index"
    tuple val(pair_id), path(reads)

    output:
    tuple val("${pair_id}"), path("${pair_id}.bam"), path("${pair_id}.junction")

    script:
    """
    STAR \
		--runThreadN ${task.cpus} \
		--twopassMode Basic \
		--readFilesCommand zcat \
		--outSAMtype BAM SortedByCoordinate \
		--outSAMstrandField intronMotif \
		--readFilesIn ${reads[0]} ${reads[1]} \
		--genomeDir STAR_genome_index \
        --chimSegmentMin 10 \
        --chimOutType Junctions \
		--outFileNamePrefix ${pair_id}
        mv ${pair_id}Aligned.sortedByCoord.out.bam ${pair_id}.bam
        mv ${pair_id}Chimeric.out.junction ${pair_id}.junction
    """
}




///////////////////////////////////////////////////////////////
// CIRCexplorer2

process CIRCEXP2_generate_annotation {
    container 'ndatth/rna-tools:v0.0.0'
    publishDir "${params.trace_dir}/circexp2", mode: 'symlink', overwrite: true
    memory '8 GB'

    input:
    path "annotation.gtf"

    output:
    path "circexp2_annotation.txt"

    script:
    """
    circexp2_generate_annotation.sh annotation.gtf circexp2_annotation.txt
    """
}


process CIRCEXP2_pipeline {
    container 'ndatth/rna-tools:v0.0.0'
    publishDir "${params.outdir}/circexp2", mode: 'copy', overwrite: true
    memory '2 GB'
    
    input:
    path "genome.fa"
    path "circexp2_annotation.txt"
    tuple val(pair_id), path("${pair_id}.bam"), path("${pair_id}.junction")

    output:
    path("${pair_id}")

    script:
    """
    CIRCexplorer2 parse -t STAR ${pair_id}.junction > CIRCexplorer2_parse.log
    CIRCexplorer2 annotate -r circexp2_annotation.txt -g genome.fa -b back_spliced_junction.bed -o ${pair_id}
    """
}


process CIRCEXP2_merge_samples { 
    container 'ndatth/rna-tools:v0.0.0'
    publishDir "${params.outdir}/circexp2", mode: 'copy', overwrite: true

    input:
    path "libsize.tsv"
    val bsj_filter
    val exp_prop
    path circexp2_outputs
 
    output:
    path "circexp2_res_merged_bsj_${bsj_filter}.tsv"
 
    script:
    """
    circexp2_merge_results.R circexp2_res_merged_bsj_${bsj_filter}.tsv libsize.tsv ${bsj_filter} ${exp_prop} ${circexp2_outputs}
    """
}



process CIRCEXP2_quantile_norm { 
    container 'ndatth/rna-tools:v0.0.0'
    publishDir "${params.outdir}/quantile_norm", mode: 'copy', overwrite: true
    memory '8 GB'
    
    input:
    path "res_merged.tsv"
    path "meta.csv"
    
    output:
    path "circexp2_quantile_norm.tsv"
 
    script:
    """
    quantile_norm.R res_merged.tsv meta.csv circexp2_quantile_norm.tsv
    """
}


process CIRCEXP2_covariate_processing { 
    
    publishDir "${params.trace_dir}/circexp2_qtl_input", mode: 'symlink', overwrite: true
    container 'ndatth/qtl-package:v0.0.0'
    memory '4 GB'

    input:
    path meta
    path exp_quatile_normalized
    path genotype_PCs 
    val n
 
    output:
    tuple val("${n}"), path("${n}_peer_covariates.tsv")
 
    script:
    """
    ## covariate_processing.R meta.csv exp_quatile_normalized.tsv genotype_PCs.tsv number_hidden_factors out_fn
    covariate_processing.R ${meta} ${exp_quatile_normalized} ${genotype_PCs} ${n} ${n}_peer_covariates.tsv
    """
}


process CIRCEXP2_chrom_splitting {

    publishDir "${params.outdir}/circexp2_qtl_input", mode: 'copy', overwrite: true
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



process CIRCEXP2_qtl_mapping {

    publishDir "${params.trace_dir}/circexp2_qtl_mapping", mode: 'symlink', overwrite: true
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


process CIRCEXP2_merge_qtl_results {

    publishDir "${params.trace_dir}/circexp2_qtl_mapping_merged", mode: 'symlink', overwrite: true
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




process CIRCEXP2_apply_qvalue {

    publishDir "${params.outdir}/qtl_mapping_apply_qvalue", mode: 'copy', overwrite: true
    container 'ndatth/qtl-package:v0.0.0'
    memory '8 GB'

    input:
    val fdr
    path merge_qtl_results

    output:
    
    tuple path("circexp2"), path("circexp2.tsv"), path("circexp2.stat"), path("circexp2.png"), path("circexp2_raw.tsv")
 
    script:
    """
    fastQTL_apply_qvalue.R circexp2 ${fdr} ${merge_qtl_results}
    """
}

process CIRCEXP2_export_all_peer_covariates { 
    
    publishDir "${params.trace_dir}/circexp2_qtl_covariates", mode: 'symlink', overwrite: true
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


process CIRCEXP2_export_optimal_peer_covariates {

    publishDir "${params.trace_dir}/circexp2_qtl_covariates", mode: 'symlink', overwrite: true
    container 'ndatth/qtl-package:v0.0.0'
    memory '2 GB'

    input:
    tuple path("circexp2"), path("circexp2.tsv"), path("circexp2.stat"), path("circexp2.png"), path("circexp2_raw.tsv")
    path all_peer_factor_covariates

    output:
    
    path("optimal_covariates.tsv") 
 
    script:
    """
    n="\$(cat circexp2)"
    cat \${n}_peer_covariates_exported.tsv > optimal_covariates.tsv
    """
}


process CIRCEXP2_qtl_mapping_nominal {
    publishDir "${params.outdir}/circexp2_qtl_mapping_nominal", mode: 'copy', overwrite: true
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


process CIRCEXP2_merge_qtl_results_nominal {

    publishDir "${params.trace_dir}/circexp2_qtl_mapping_nominal", mode: 'symlink', overwrite: true
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





///////////////////////////////////////////////////////////////
// Salmon

process SALMON_index_cdna {
    container 'ndatth/rna-tools:v0.0.0'

    input:
    path "cdna.fa"

    output:
    path "index_cdna"

    script:
    """
    salmon index -t cdna.fa -i index_cdna
    """
}

process SALMON_quantify {
    container 'ndatth/rna-tools:v0.0.0'
    publishDir "${params.outdir}/salmon", mode: 'copy', overwrite: true
    cpus 32
    memory '32 GB'

    input:
    path "index_cdna"
    tuple val(pair_id), path(reads)

    output:
    path "${pair_id}"

    script:
    """
    salmon quant -i index_cdna -l A -1 ${reads[0]} -2 ${reads[1]} -p ${task.cpus} --validateMappings -o ${pair_id}_quant
    mv ${pair_id}_quant/quant.sf ${pair_id}
    """
}


process SALMON_merge_samples {
    container 'ndatth/rna-tools:v0.0.0'
    publishDir "${params.outdir}/salmon", mode: 'copy', overwrite: true
    memory '8 GB'

    input:
    path "annotation.gtf"
    path salmon_outputs

    output:
    path "salmon_res_merged.tsv"

    script:
    """
    salmon_merge_results.R annotation.gtf salmon_res_merged.tsv ${salmon_outputs}
    """
}


process SALMON_quantile_norm { 
    container 'ndatth/rna-tools:v0.0.0'
    publishDir "${params.outdir}/quantile_norm", mode: 'copy', overwrite: true
    memory '8 GB'
    
    input:
    path "res_merged.tsv"
    path "meta.csv"
    val exp_prop
    path "annotation.gtf"
 
    output:
    path "salmon_quantile_norm.tsv"
 
    script:
    """
    extract_gene_gtf.py annotation.gtf gene_info.tsv
    quantile_norm_salmon.R res_merged.tsv meta.csv gene_info.tsv ${exp_prop} salmon_quantile_norm.tsv
    """
}


process SALMON_covariate_processing { 
    
    publishDir "${params.trace_dir}/salmon_qtl_input"
    container 'ndatth/qtl-package:v0.0.0'
    memory '4 GB'

    input:
    path meta
    path exp_quatile_normalized
    path genotype_PCs 
    val n
 
    output:
    tuple val("${n}"), path("${n}_peer_covariates.tsv")
 
    script:
    """
    ## covariate_processing.R meta.csv exp_quatile_normalized.tsv genotype_PCs.tsv number_hidden_factors out_fn
    covariate_processing.R ${meta} ${exp_quatile_normalized} ${genotype_PCs} ${n} ${n}_peer_covariates.tsv
    """
}


process SALMON_chrom_splitting {

    publishDir "${params.outdir}/salmon_qtl_input", mode: 'copy', overwrite: true
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



process SALMON_qtl_mapping {

    publishDir "${params.trace_dir}/salmon_qtl_mapping", mode: 'symlink', overwrite: true
    container 'ndatth/qtl-package:v0.0.0'
    cpus 32
    memory '45 GB'

    input:
    tuple val("n"), path("covariates.tsv"), val("chr"), path("chr.vcf.gz"), path("chr.vcf.gz.tbi"), path("chr.bed.gz"), path("chr.bed.gz.tbi")
    val window_size

    output:
    
    tuple val("${n}"), path("peer_${n}_chr_${chr}_permute_pass.tsv.gz")
 
    script:
    """
    #!/bin/bash

    for i in \$(seq 1 ${task.cpus}); do
        fastQTL --vcf chr.vcf.gz --bed chr.bed.gz --cov covariates.tsv --window ${window_size} --permute 1000 10000 --chunk \${i} ${task.cpus} --out peer_${n}_chr_${chr}chunk_\${i}_permute_pass.tsv &
    done
    wait
    cat peer_${n}_chr_${chr}chunk_*_permute_pass.tsv | gzip -c > peer_${n}_chr_${chr}_permute_pass.tsv.gz
    """
}



process SALMON_merge_qtl_results {

    publishDir "${params.trace_dir}/salmon_qtl_mapping_merged", mode: 'symlink', overwrite: true
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


process SALMON_apply_qvalue {

    publishDir "${params.outdir}/qtl_mapping_apply_qvalue", mode: 'copy', overwrite: true
    container 'ndatth/qtl-package:v0.0.0'
    memory '8 GB'

    input:
    val fdr
    path merge_qtl_results

    output:
    
    tuple path("salmon"), path("salmon.tsv"), path("salmon.stat"), path("salmon.png"), path("salmon_raw.tsv")
 
    script:
    """
    fastQTL_apply_qvalue.R salmon ${fdr} ${merge_qtl_results}
    """
}

process SALMON_export_all_peer_covariates { 
    
    publishDir "${params.trace_dir}/salmon_qtl_covariates", mode: 'symlink', overwrite: true
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


process SALMON_export_optimal_peer_covariates {

    publishDir "${params.trace_dir}/salmon_qtl_covariates", mode: 'symlink', overwrite: true
    container 'ndatth/qtl-package:v0.0.0'
    memory '2 GB'

    input:
    tuple path("salmon"), path("salmon.tsv"), path("salmon.stat"), path("salmon.png") , path("salmon_raw.tsv")
    path all_peer_factor_covariates

    output:
    
    path("optimal_covariates.tsv") 
 
    script:
    """
    n="\$(cat salmon)"
    cat \${n}_peer_covariates_exported.tsv > optimal_covariates.tsv
    """
}


process SALMON_qtl_mapping_nominal {

    publishDir "${params.outdir}/salmon_qtl_mapping_nominal", mode: 'copy', overwrite: true
    container 'ndatth/qtl-package:v0.0.0'
    memory '4 GB'

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


process SALMON_merge_qtl_results_nominal {

    publishDir "${params.trace_dir}/salmon_qtl_mapping_nominal", mode: 'symlink', overwrite: true
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




