

process EXTRACT_SEQVARS {

    publishDir "${params.outputs_dir}/${params.run_id}/${vtype.toLowerCase()}", mode: 'symlink'

    input:
    tuple val(sample), val(meta), path(vcf), path(scna)
    val vtype

    output:
    path "${sample}.${vtype.toLowerCase()}.tsv"

    script:
    def allow_utr       = params.seqvars.allow_utr       ? '--allow-utr'       : ''
    def allow_splice    = params.seqvars.allow_splice    ? '--allow-splice'    : ''
    def allow_noncoding = params.seqvars.allow_noncoding ? '--allow-noncoding' : ''
    """
    python ${params.scripts.extract_seqvars} \
        --vcf ${vcf} \
        --scna ${scna} \
        --purity ${meta.purity} \
        --vtype ${vtype} \
        ${allow_utr} \
        ${allow_splice} \
        ${allow_noncoding} \
        --allow-intergenic-dist ${params.seqvars.allow_intergenic_dist} \
        --outfile ${sample}.${vtype.toLowerCase()}.tsv

    """

}


process EXTRACT_STRUCTVARS {

    publishDir "${params.outputs_dir}/${params.run_id}/sv", mode: 'symlink'

    input:
    tuple val(sample), val(meta), path(vcf)

    output:
    path "${sample}.sv.tsv"

    script:
    """
    python ${params.scripts.extract_structvars} \
        --vcf ${vcf} \
        --outfile ${sample}.sv.tsv
    """

}

process EXTRACT_CNA {

    publishDir "${params.outputs_dir}/${params.run_id}/cna", mode: 'symlink'

    input:
    tuple val(sample), val(meta), path(scna)
    path gff

    output:
    path "${sample}.cna.tsv"

    script:
    def wgd               = meta.wgd                         ? '--wgd' : ''
    def allow_subclonal   = params.cnavars.allow_subclonal   ? '--allow-subclonal' : ''
    def allow_amp         = params.cnavars.allow_amp         ? '--allow-amp' : ''
    def allow_deep_del    = params.cnavars.allow_deep_del    ? '--allow-deep-del' : ''
    def allow_shallow_del = params.cnavars.allow_shallow_del ? '--allow-shallow-del' : ''
    def allow_noncoding   = params.cnavars.allow_noncoding   ? '--allow-noncoding' : ''
    def allow_chrX        = params.cnavars.allow_chrX        ? '--allow-x' : ''
    """
    python ${params.scripts.extract_cna} \
        --scna ${scna} \
        --gff ${gff} \
        --min-span ${params.cnavars.min_span} \
        ${wgd} \
        ${allow_subclonal} \
        ${allow_amp} \
        ${allow_deep_del} \
        ${allow_shallow_del} \
        ${allow_noncoding} \
        ${allow_chrX} \
        --outfile ${sample}.cna.tsv
    """

}


process MERGE_DONORS {

    publishDir "${params.outputs_dir}/${params.run_id}/merged_mutations", mode: 'copy'

    input:
    path svfiles, stageAs: 'svs/*'
    path cnafiles, stageAs: 'cna/*'
    path snvfiles, stageAs: 'snvs/*'
    path indelfiles, stageAs: 'indels/*'

    output:
    path "merged.tsv", emit: merged
    path "merged.log", emit: log

    script:
    """
    python ${params.scripts.merge_donors} \
        --svdir svs \
        --cnadir cna \
        --snvdir snvs \
        --indeldir indels \
        --outfile merged.tsv \
        --logfile merged.log
    """

}

