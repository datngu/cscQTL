#!/usr/bin/env nextflow

// include { RECOUNT_get_consensus; } from './module/recount'
// RECOUNT_covariate_processing; RECOUNT_chrom_splitting; RECOUNT_qtl_mapping; RECOUNT_merge_qtl_results; RECOUNT_apply_qvalue; RECOUNT_export_all_peer_covariates;


nextflow.enable.dsl=2

process RENAME_circall {
    container 'ndatth/qtl-package:v0.0.0'
    publishDir "${params.trace_dir}/rename_circ_pipelines", mode: 'symlink', overwrite: true
    memory '8 GB'
    cpus 1

    input:
    path ID

    output:
    path "${ID}.circall"

    script:
    """
    cat ${ID} > "${ID}.circall"
    """
}


process RENAME_ciri2 {
    container 'ndatth/qtl-package:v0.0.0'
    publishDir "${params.trace_dir}/rename_circ_pipelines", mode: 'symlink', overwrite: true
    memory '8 GB'
    cpus 1

    input:
    path ID

    output:
    path "${ID}.ciri2"

    script:
    """
    cat ${ID} > "${ID}.ciri2"
    """
}


process RENAME_circexp2 {
    container 'ndatth/qtl-package:v0.0.0'
    publishDir "${params.trace_dir}/rename_circ_pipelines", mode: 'symlink', overwrite: true
    memory '8 GB'
    cpus 1

    input:
    path ID

    output:
    path "${ID}.circexp2"

    script:
    """
    cat ${ID} > "${ID}.circexp2"
    """
}


process RECOUNT_get_consensus {
    container 'ndatth/qtl-package:v0.0.0'
    publishDir "${params.outdir}/recount_pseudo", mode: 'copy', overwrite: true
    memory '8 GB'
    cpus 1

    input:
    path gtf
    path circall_res
    path ciri2_res
    path circexp2_res

    output:
    path "circRNA_consensus.txt"

    script:
    """
    get_circRNA_consensus.py --gtf $gtf --circall $circall_res --ciri2 $ciri2_res --circexp $circexp2_res --consensus $params.consensus --out circRNA_consensus.txt --cutoff $params.bsj_filter
    """
}


process RECOUNT_get_peseudo_seqs {
    container 'ndatth/rna-tools:v0.0.0'
    publishDir "${params.outdir}/recount_pseudo", mode: 'copy', overwrite: true
    memory '8 GB'
    cpus 1

    input:
    path genome
    path circRNA

    output:
    path "circRNA_consensus.fa"

    script:
    """
    generate_peseudo_BSJ_sequences.R genome=$genome circRNA=$circRNA out=circRNA_consensus.fa
    """
}


process RECOUNT_index_peseudo_seqs {
    container 'ndatth/rna-tools:v0.0.0'
    publishDir "${params.trace_dir}/recount_pseudo_idx", mode: 'symlink', overwrite: true
    memory '8 GB'
    cpus 1

    input:
    path circRNA_seq

    output:
    path "index_pseudo_seqs"

    script:
    """
    TxIndexer -t $circRNA_seq -o index_pseudo_seqs
    """
}


process RECOUNT_mapping {
    container 'ndatth/rna-tools:v0.0.0'
    publishDir "${params.outdir}/recount_mapping", mode: 'copy', overwrite: true
    memory '32 GB'
    cpus 32

    input:
    path "index_cdna"
    path "index_bsj" 
    tuple val(pair_id), path(reads)

    output:
    tuple val("${pair_id}"), path("${pair_id}_Circall_wt_fragmentInfo.txt"), path("${pair_id}_BSJ_mapped_read.txt")

    script:
    """
    Circall_recount.sh -txIdx index_cdna -bsjIdx index_bsj -read1 ${reads[0]} -read2 ${reads[1]} -o ${pair_id} -p ${task.cpus}
    """
}


process RECOUNT_count {
    container 'ndatth/rna-tools:v0.0.0'
    publishDir "${params.outdir}/recount", mode: 'copy', overwrite: true
    memory '4 GB'

    input:
    tuple val(pair_id), path(fragmentInfo), path(bsj_mapped_read)

    output:
    path "${pair_id}.recount"

    script:
    """
    recount_BSJ.py --bsj "${bsj_mapped_read}" --fragment --hit 50 --anchor 7 --info "${fragmentInfo}" --out "${pair_id}.recount"
    """
}



process RECOUNT_merge_samples {
    container 'ndatth/rna-tools:v0.0.0'
    publishDir "${params.outdir}/recount", mode: 'copy', overwrite: true
    memory '8 GB'

    input:
    val bsj_filter
    val exp_prop
    path recount_outputs

    output:
    path "recount_res_merged_bsj_${bsj_filter}.tsv"

    script:
    """
    recount_merge_results.R "recount_res_merged_bsj_${bsj_filter}.tsv" ${bsj_filter} ${exp_prop} ${recount_outputs}
    """
}


process RECOUNT_quantile_norm { 
    container 'ndatth/rna-tools:v0.0.0'
    publishDir "${params.outdir}/quantile_norm", mode: 'copy', overwrite: true
    memory '8 GB'
    
    input:
    path "res_merged.tsv"
    path "meta.csv"
    
    output:
    path "recount_quantile_norm.tsv"
 
    script:
    """
    quantile_norm.R res_merged.tsv meta.csv recount_quantile_norm.tsv
    """
}


process RECOUNT_covariate_processing { 
    
    publishDir "${params.trace_dir}/recount_qtl_input", mode: 'symlink', overwrite: true
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



process RECOUNT_chrom_splitting {

    publishDir "${params.outdir}/recount_qtl_input", mode: 'copy', overwrite: true
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



process RECOUNT_qtl_mapping {

    publishDir "${params.trace_dir}/recount_qtl_mapping", mode: 'symlink', overwrite: true
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



process RECOUNT_merge_qtl_results {

    publishDir "${params.trace_dir}/recount_qtl_mapping_merged", mode: 'symlink', overwrite: true
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


process RECOUNT_apply_qvalue {

    publishDir "${params.outdir}/qtl_mapping_apply_qvalue", mode: 'copy', overwrite: true
    container 'ndatth/qtl-package:v0.0.0'
    memory '8 GB'

    input:
    val fdr
    path merge_qtl_results

    output:
    
    tuple path("recount"), path("recount.tsv"), path("recount.stat"), path("recount.png"), path("recount_raw.tsv") 
 
    script:
    """
    fastQTL_apply_qvalue.R recount ${fdr} ${merge_qtl_results}
    """
}


process RECOUNT_export_all_peer_covariates { 
    
    publishDir "${params.trace_dir}/recount_qtl_covariates", mode: 'symlink', overwrite: true
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



process RECOUNT_export_optimal_peer_covariates {

    publishDir "${params.trace_dir}/recount_qtl_covariates", mode: 'symlink', overwrite: true
    container 'ndatth/qtl-package:v0.0.0'
    memory '2 GB'

    input:
    tuple path("recount"), path("recount.tsv"), path("recount.stat"), path("recount.png"), path("recount_raw.tsv") 
    path all_peer_factor_covariates

    output:
    
    path("optimal_covariates.tsv") 
 
    script:
    """
    n="\$(cat recount)"
    cat \${n}_peer_covariates_exported.tsv > optimal_covariates.tsv
    """
}



process RECOUNT_qtl_mapping_nominal {

    publishDir "${params.outdir}/recount_qtl_mapping_nominal", mode: 'copy', overwrite: true
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




process RECOUNT_merge_qtl_results_nominal {

    publishDir "${params.trace_dir}/recount_qtl_mapping_nominal", mode: 'symlink', overwrite: true
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

