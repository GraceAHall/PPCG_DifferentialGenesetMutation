

process SNPEFF_DOWNLOAD {

    debug true
    cpus 1
    time "20m"
    memory "4 GB"
    publishDir "${params.outputs_dir}/${params.run_id}/snpeff_download", mode: 'copy'
    container 'quay.io/biocontainers/snpeff:5.3.0a--hdfd78af_1'
    
    output:
    path "snpeff_data/GRCh37.75"
    
    script:
    """
    mkdir snpeff_data
    snpEff download -v GRCh37.75
    mv /usr/local/share/snpeff-5.3.0a-1/./data/GRCh37.75 snpeff_data
    """
}

process SNPEFF_RUN {

    debug true
    cpus 1
    time "20m"
    memory "8 GB"
    publishDir "${params.outputs_dir}/${params.run_id}/snpeff_run", mode: 'copy'
    container 'quay.io/biocontainers/snpeff:5.3.0a--hdfd78af_1'
    
    input:
    tuple val(sample), val(meta), path(vcf)
    path data_dir
    
    output:
    tuple val(sample), val(meta), path("${sample}.snpEff.vcf")
    
    script:
    """
    sed 's/chrM/chrMT/' ${vcf} > tmp.vcf
    snpEff ann -canon -dataDir ${data_dir} GRCh37.75 tmp.vcf > ${sample}.snpEff.vcf
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


process MERGE_VARIANTS {

    publishDir "${params.outputs_dir}/${params.run_id}/merged_variants", mode: 'copy'

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
    python ${params.scripts.merge_variants} \
        --svdir svs \
        --cnadir cna \
        --snvdir snvs \
        --indeldir indels \
        --outfile merged.tsv \
        --logfile merged.log
    """

}


// does the CADD feature annotation mirror the VCF ANNOVAR feature annotation? No.
// ANNOVAR and CADD use different transcripts. 
// CADD scores should be used with caution. 
// Eg. PPCG0086a - chr1:9664540 G>T (TMEM201). 
//   - ANNOVAR says 3_prime_UTR_variant, CADD says intron_variant
process FILTER_VARIANTS_COMBIMETS {

    publishDir "${params.outputs_dir}/${params.run_id}/filtered_variants", mode: 'copy'

    input:
    path merged_variants
    path trees_dir
    path ccfs_dir
    path mettraj_clones
    path samplesheet

    output:
    path "filtered.tsv"

    script:
    """
    python ${params.scripts.filter_variants_combimets} \
        --mutations ${merged_variants} \
        --trees-dir ${trees_dir} \
        --ccfs-dir ${ccfs_dir} \
        --mettraj-clones ${mettraj_clones} \
        --samplesheet ${samplesheet} \
        --zscore-thresh ${params.filter_variants.hypermutator_zscore} \
        --outfile-mutations filtered.tsv \
        --outfile-summary hypermutators.tsv
    """

}

process FILTER_VARIANTS_GENERIC {

    publishDir "${params.outputs_dir}/${params.run_id}/filtered_variants", mode: 'copy'

    input:
    path merged_variants

    output:
    path "filtered.tsv"

    script:
    """
    python ${params.scripts.filter_variants_generic} \
        --mutations ${merged_variants} \
        --zscore-thresh ${params.filter_variants.hypermutator_zscore} \
        --outfile-mutations filtered.tsv \
        --outfile-summary hypermutators.tsv
    """

}

