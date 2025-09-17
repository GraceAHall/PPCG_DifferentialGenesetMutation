

process LOG_SETTINGS {
    
    publishDir "${params.outputs_dir}/${params.run_id}/pipeline_info", mode: 'copy'
    
    output:
    path "settings.txt"
    
    script:
    """
    echo "Nextflow Config Settings:" > settings.txt
    echo "=========================" >> settings.txt
    echo "Nextflow version: ${nextflow.version}" >> settings.txt
    echo "Launch directory: ${workflow.launchDir}" >> settings.txt
    echo "Work directory: ${workflow.workDir}" >> settings.txt
    echo "" >> settings.txt
    
    echo "Parameters:" >> settings.txt
    echo "===========" >> settings.txt

    echo "SNV|INDEL extraction" >> settings.txt
    echo "- allow_utr: ${params.seqvars.allow_utr}" >> settings.txt
    echo "- allow_splice: ${params.seqvars.allow_splice}" >> settings.txt
    echo "- allow_noncoding: ${params.seqvars.allow_noncoding}" >> settings.txt
    echo "- allow_intergenic_dist: ${params.seqvars.allow_intergenic_dist}" >> settings.txt
    echo "" >> settings.txt
    
    echo "CNA extraction" >> settings.txt
    echo "- min_span: ${params.cnavars.min_span}" >> settings.txt
    echo "- allow_amp: ${params.cnavars.allow_amp}" >> settings.txt
    echo "- allow_deep_del: ${params.cnavars.allow_deep_del}" >> settings.txt
    echo "- allow_shallow_del: ${params.cnavars.allow_shallow_del}" >> settings.txt
    echo "- allow_noncoding: ${params.cnavars.allow_noncoding}" >> settings.txt
    echo "- allow_chrX: ${params.cnavars.allow_chrX}" >> settings.txt
    echo "" >> settings.txt
    
    echo "Geneset selection" >> settings.txt
    echo "- min_genes: ${params.select_genesets.min_genes}" >> settings.txt
    echo "- max_genes: ${params.select_genesets.max_genes}" >> settings.txt
    echo "- min_genes_hitprop: ${params.select_genesets.min_genes_hitprop}" >> settings.txt
    echo "- min_cohort_hitprop: ${params.select_genesets.min_cohort_hitprop}" >> settings.txt
    echo "- max_gene_reliance: ${params.select_genesets.max_gene_reliance}" >> settings.txt
    echo "" >> settings.txt

    echo "Background enrichment" >> settings.txt
    echo "- budget: ${params.background_enrichment.budget}" >> settings.txt
    echo "" >> settings.txt
    
    echo "Results summary" >> settings.txt
    echo "- max_plots: ${params.summarise_results.max_plots}" >> settings.txt
    echo "" >> settings.txt
    """
}
