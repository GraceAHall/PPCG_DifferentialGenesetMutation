

process SNPEFF_DOWNLOAD {

    debug true
    cpus 1
    time "20m"
    memory "4 GB"
    publishDir "${params.outputs_dir}/${params.run_id}/variant_extraction/snpeff_download", mode: 'copy'
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
    publishDir "${params.outputs_dir}/${params.run_id}/variant_extraction/snpeff_run", mode: 'copy'
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

    publishDir "${params.outputs_dir}/${params.run_id}/variant_extraction/sv", mode: 'symlink'

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

    publishDir "${params.outputs_dir}/${params.run_id}/variant_extraction/${vtype.toLowerCase()}", mode: 'symlink'

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

    publishDir "${params.outputs_dir}/${params.run_id}/variant_extraction/cna", mode: 'symlink'

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

    publishDir "${params.outputs_dir}/${params.run_id}/variant_processing", mode: 'copy'

    input:
    path svfiles, stageAs: 'svs/*'
    path cnafiles, stageAs: 'cna/*'
    path snvfiles, stageAs: 'snvs/*'
    path indelfiles, stageAs: 'indels/*'

    output:
    path "mutations.merged.tsv", emit: mutations
    path "mutations.merged.log", emit: log

    script:
    """
    python ${params.scripts.merge_variants} \
        --svdir svs \
        --cnadir cna \
        --snvdir snvs \
        --indeldir indels \
        --outfile mutations.merged.tsv \
        --logfile mutations.merged.log
    """

}

process STANDARDISE_AND_FILTER_VARIANTS {

    publishDir "${params.outputs_dir}/${params.run_id}/variant_processing", mode: 'symlink'

    input:
    path mutations
    path genesets
    path gff
    path hgnc

    output:
    path 'mutations.filtered.tsv', emit: 'mutations'
    path 'genesets.tsv', emit: 'genesets'
    path 'sizes.tsv', emit: 'sizes'

    script:
    """
    python ${params.scripts.standardise_filter_variants} \
        --mutations ${mutations} \
        --genesets ${genesets} \
        --gff ${gff} \
        --hgnc ${hgnc} \
        --outfile-muts mutations.filtered.tsv \
        --outfile-gsets genesets.tsv \
        --outfile-sizes sizes.tsv

    """

}

// does the CADD feature annotation mirror the VCF ANNOVAR feature annotation? No.
// ANNOVAR and CADD use different transcripts. 
// CADD scores should be used with caution. 
// Eg. PPCG0086a - chr1:9664540 G>T (TMEM201). 
//   - ANNOVAR says 3_prime_UTR_variant, CADD says intron_variant
process ASSIGN_CLONES {

    publishDir "${params.outputs_dir}/${params.run_id}/variant_processing", mode: 'copy'

    input:
    path merged_variants
    path trees_dir
    path ccfs_dir
    path mettraj_clones
    path samplesheet

    output:
    path "mutations.assigned.tsv"

    script:
    """
    python ${params.scripts.assign_clones} \
        --mutations ${merged_variants} \
        --trees-dir ${trees_dir} \
        --ccfs-dir ${ccfs_dir} \
        --mettraj-clones ${mettraj_clones} \
        --samplesheet ${samplesheet} \
        --outfile-mutations mutations.assigned.tsv
    """

}

process ASSIGN_CLONES_PLACEHOLDER {

    publishDir "${params.outputs_dir}/${params.run_id}/variant_processing", mode: 'copy'

    input:
    path merged_variants

    output:
    path "mutations.assigned.tsv"

    script:
    """
    python ${params.scripts.assign_clones_placeholder} \
        --mutations ${merged_variants} \
        --outfile-mutations mutations.assigned.tsv
    """

}

