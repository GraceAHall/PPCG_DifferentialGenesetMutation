
process HARMONISE_DATA {

    publishDir "${params.outputs_dir}/${params.run_id}/harmonisation", mode: 'symlink'

    input:
    path mutations
    path genesets
    path sizes 

    output:
    tuple path(mutations), path('genesets.tsv'), path('sizes.tsv')

    script:
    """
    python ${params.scripts.harmonise_data} \
        --mutations ${mutations} \
        --genesets ${genesets} \
        --sizes ${sizes} \
        --gff ${params.refseq_gff_path} \
        --gtf ${params.ensembl_gtf_path} \
        --hgnc ${params.hgnc_path} \
        --outfile-sizes sizes.tsv \
        --outfile-gsets genesets.tsv

    """

}

process SPLIT_MUTATIONS {
    
    publishDir "${params.outputs_dir}/${params.run_id}/${posclass}/split_mutations", mode: 'symlink'
    
    input:
    tuple val(posclass), path(mutations), path(genesets), path(sizes), path(samplesheet)

    output:
    tuple val(posclass), path('posmuts.tsv'), path('negmuts.tsv'), path(genesets), path(sizes)

    script:
    """
    python ${params.scripts.split_mutations} \
        --posclass ${posclass} \
        --mutations ${mutations} \
        --samplesheet ${samplesheet} \
        --outfile-pos posmuts.tsv \
        --outfile-neg negmuts.tsv

    """

}

process SELECT_GENESETS {
    
    publishDir "${params.outputs_dir}/${params.run_id}/${posclass}/select_genesets", mode: 'symlink'
    
    input:
    tuple val(posclass), path(posmuts), path(negmuts), path(genesets), path(sizes)

    output:
    tuple val(posclass), path('selection.tsv')

    script:
    """
    python ${params.scripts.select_genesets} \
        --posmuts ${posmuts} \
        --negmuts ${negmuts} \
        --genesets ${genesets} \
        --min-genes ${params.select_genesets.min_genes} \
        --max-genes ${params.select_genesets.max_genes} \
        --min-genes-hitprop ${params.select_genesets.min_genes_hitprop} \
        --min-cohort-hitprop ${params.select_genesets.min_cohort_hitprop} \
        --max-gene-reliance ${params.select_genesets.max_gene_reliance} \
        --outfile selection.tsv
    """

}

process DIFFERENTIAL_MUTATION {
    
    publishDir "${params.outputs_dir}/${params.run_id}/${posclass}/differential_mutation", mode: 'symlink'
    
    input:
    tuple val(posclass), path(posmuts), path(negmuts), path(genesets), path(sizes), path(selection)

    output:
    tuple val(posclass), path('diffmut.tsv')

    script:
    """
    python ${params.scripts.differential_mutation} \
        --posmuts ${posmuts} \
        --negmuts ${negmuts} \
        --genesets ${genesets} \
        --selection ${selection} \
        --outfile diffmut.tsv
    """

}

process BACKGROUND_ENRICHMENT {
    
    publishDir "${params.outputs_dir}/${params.run_id}/${posclass}/background_enrichment", mode: 'symlink'
    
    input:
    tuple val(posclass), path(posmuts), path(negmuts), path(genesets), path(sizes), path(selection)

    output:
    tuple val(posclass), path('pos.enrich.tsv'), path('neg.enrich.tsv')

    script:
    """
    python ${params.scripts.background_enrichment} \
        --mutations ${posmuts} \
        --genesets ${genesets} \
        --selection ${selection} \
        --sizes ${sizes} \
        --budget ${params.background_enrichment.budget} \
        --outfile pos.enrich.tsv
    
    python ${params.scripts.background_enrichment} \
        --mutations ${negmuts} \
        --genesets ${genesets} \
        --selection ${selection} \
        --sizes ${sizes} \
        --budget ${params.background_enrichment.budget} \
        --outfile neg.enrich.tsv
    """

}

process SUMMARISE_RESULTS {
    
    publishDir "${params.outputs_dir}/${params.run_id}/${posclass}/results", mode: 'copy'
    
    input:
    tuple val(posclass), path(posmuts), path(negmuts), path(genesets), path(sizes), path(selection), path(diffmut), path(posenrich), path(negenrich)

    output:
    path "${posclass}.tsv"
    path "*.png"

    script:
    """
    python ${params.scripts.summarise_results} \
        --posclass ${posclass} \
        --genesets ${genesets} \
        --selection ${selection} \
        --posmuts ${posmuts} \
        --negmuts ${negmuts} \
        --diffmut ${diffmut} \
        --posenrich ${posenrich} \
        --negenrich ${negenrich} \
        --max-plots ${params.summarise_results.max_plots}
    """

}

