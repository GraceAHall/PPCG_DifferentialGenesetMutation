
// process HARMONISE_DATA {

//     publishDir "${params.outputs_dir}/${params.run_id}/contrasts/harmonisation", mode: 'symlink'

//     input:
//     path mutations
//     path genesets
//     path sizes 

//     output:
//     tuple path(mutations), path('genesets.tsv'), path('sizes.tsv')

//     script:
//     """
//     python ${params.scripts.harmonise_data} \
//         --mutations ${mutations} \
//         --genesets ${genesets} \
//         --sizes ${sizes} \
//         --gff ${params.refseq_gff_path} \
//         --gtf ${params.ensembl_gtf_path} \
//         --hgnc ${params.hgnc_path} \
//         --outfile-sizes sizes.tsv \
//         --outfile-gsets genesets.tsv

//     """

// }

// process PREPARE_DOWNSTREAM {

//     publishDir "${params.outputs_dir}/${params.run_id}/contrasts/prepare_downstream", mode: 'symlink'

//     input:
//     path mutations
//     path genesets
//     path gff
//     path hgnc

//     output:
//     tuple path('mutations.tsv'), path('genesets.tsv'), path('sizes.tsv')

//     script:
//     """
//     python ${params.scripts.prepare_downstream} \
//         --mutations ${mutations} \
//         --genesets ${genesets} \
//         --gff ${gff} \
//         --hgnc ${hgnc} \
//         --zscore-thresh ${params.filter_variants.hypermutator_zscore} \
//         --outfile-muts mutations.tsv \
//         --outfile-gsets genesets.tsv \
//         --outfile-sizes sizes.tsv

//     """

// }

process SPLIT_MUTATIONS {
    
    publishDir "${params.outputs_dir}/${params.run_id}/contrasts/${posclass}/split_mutations", mode: 'symlink'
    
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
    
    publishDir "${params.outputs_dir}/${params.run_id}/contrasts/${posclass}/select_genesets", mode: 'symlink'
    
    input:
    tuple val(posclass), path(posmuts), path(negmuts), path(genesets), path(sizes)

    output:
    tuple val(posclass), path('selection.tsv')

    script:
    """
    echo "hi"
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


process SINGLE_GENE_ANALYSIS {
    
    publishDir "${params.outputs_dir}/${params.run_id}/contrasts/${posclass}/single_gene_analysis", mode: 'symlink'
    
    input:
    tuple val(posclass), path(posmuts), path(negmuts), path(genesets), path(sizes)

    output:
    tuple path('summary.log'), path('summary.tsv')

    script:
    """
    echo "hi"
    python ${params.scripts.single_gene_analysis} \
        --posmuts ${posmuts} \
        --negmuts ${negmuts} \
        --outfile summary
    """

}

process DIFFERENTIAL_MUTATION {
    
    publishDir "${params.outputs_dir}/${params.run_id}/contrasts/${posclass}/differential_mutation", mode: 'symlink'
    
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
    
    publishDir "${params.outputs_dir}/${params.run_id}/contrasts/${posclass}/background_enrichment", mode: 'symlink'
    
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
        --outfile pos.enrich.tsv
    
    python ${params.scripts.background_enrichment} \
        --mutations ${negmuts} \
        --genesets ${genesets} \
        --selection ${selection} \
        --sizes ${sizes} \
        --outfile neg.enrich.tsv
    """

}

process SUMMARISE_RESULTS {
    
    publishDir "${params.outputs_dir}/${params.run_id}/contrasts/${posclass}/results", mode: 'copy'
    
    input:
    tuple val(posclass), path(posmuts), path(negmuts), path(genesets), path(sizes), path(selection), path(diffmut), path(posenrich), path(negenrich)

    output:
    path "${posclass}.full.tsv"
    // path "${posclass}.sig.tsv"
    path "oncoprints/"

    script:
    """
    echo "hh"
    mkdir oncoprints
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

