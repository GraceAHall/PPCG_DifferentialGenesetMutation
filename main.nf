
include { LOG_SETTINGS } from './modules/misc'

include { SNPEFF_DOWNLOAD } from './modules/extract'
include { SNPEFF_RUN } from './modules/extract'
include { EXTRACT_STRUCTVARS } from './modules/extract'
include { EXTRACT_SEQVARS as EXTRACT_SNVS } from './modules/extract'
include { EXTRACT_SEQVARS as EXTRACT_INDELS } from './modules/extract'
include { EXTRACT_CNA } from './modules/extract'
include { MERGE_VARIANTS } from './modules/extract'
include { STANDARDISE_AND_FILTER_VARIANTS } from './modules/extract'
include { ASSIGN_CLONES } from './modules/extract'
include { ASSIGN_CLONES_PLACEHOLDER } from './modules/extract'

// include { HARMONISE_DATA } from './modules/contrast'
// include { PREPARE_DOWNSTREAM } from './modules/contrast'
include { SPLIT_MUTATIONS } from './modules/contrast'
include { SINGLE_GENE_ANALYSIS } from './modules/contrast'
include { SELECT_GENESETS } from './modules/contrast'
include { DIFFERENTIAL_MUTATION } from './modules/contrast'
include { BACKGROUND_ENRICHMENT } from './modules/contrast'
include { SUMMARISE_RESULTS } from './modules/contrast'

// TODO: Swap CNA extraction to use GENCODE gff. 

genesets        = file(params.genesets)
hgnc            = file(params.hgnc_path)
gff             = file(params.gencode_gff_path)
cosmic          = file(params.cosmic_gene_locations)

Channel.fromPath(params.samplesheet)
    | splitCsv(header: true, sep: '\t')
    | map { row ->
        def sample = row.sample
        def metadata = [
            sample: row.sample,
            donor: row.donor,
            purity: row.purity,
            ploidy: row.ploidy, 
            wgd: row.WGD,
            cohort: row.cohort
        ]
        return tuple(sample, metadata)
    }
    | set { samples_ch }


workflow {

    LOG_SETTINGS()
    
    // --- INITAL FILE WRANGLING --- //
    ch_samplesheet  = Channel.fromPath(params.samplesheet) 
    sv_files        = Channel.fromPath("${params.basedir}/data/sv_annot/*.vcf")
    snv_files       = Channel.fromPath("${params.basedir}/data/snv/*.vcf")
    indel_files     = Channel.fromPath("${params.basedir}/data/indel/*.vcf")
    scna_files      = Channel.fromPath("${params.basedir}/data/scna/*.txt")

    all_svs     = sv_files      | map { filepath -> tuple(filepath.simpleName.toString()[0..8], filepath) }
    all_snvs    = snv_files     | map { filepath -> tuple(filepath.simpleName.toString()[0..8], filepath) }
    all_indels  = indel_files   | map { filepath -> tuple(filepath.simpleName.toString()[0..8], filepath) }
    all_scna    = scna_files    | map { filepath -> tuple(filepath.simpleName.toString()[0..8], filepath) }

    ch_indels_raw   = samples_ch | join(all_indels) | join(all_scna)
    ch_snvs_raw     = samples_ch | join(all_snvs) | join(all_scna)
    ch_svs_raw      = samples_ch | join(all_svs)
    ch_cna_raw      = samples_ch | join(all_scna)

    // --- EXTRACT VARIANTS --- //
    // extract snvs, indels, svs and cna into standardised .tsv files
    if (params.snpeff_completed) {
        ch_svs_snpeff = ch_svs_raw
    } else {
        ch_snpeff_dir = SNPEFF_DOWNLOAD()
        ch_svs_snpeff = SNPEFF_RUN(ch_svs_raw, ch_snpeff_dir)
    }
    ch_svs      = EXTRACT_STRUCTVARS(ch_svs_snpeff)
    ch_indels   = EXTRACT_INDELS(ch_indels_raw, 'INDEL')
    ch_snvs     = EXTRACT_SNVS(ch_snvs_raw, 'SNV')
    ch_cna      = EXTRACT_CNA(ch_cna_raw, cosmic)
    
    // merge variants into single table
    MERGE_VARIANTS(
        ch_svs.collect(),
        ch_cna.collect(),
        ch_snvs.collect(),
        ch_indels.collect()
    )
    ch_merged = MERGE_VARIANTS.out.mutations
    
    // filter variants & harmonise gene symbols between various data sources
    STANDARDISE_AND_FILTER_VARIANTS(ch_merged, genesets, gff, hgnc) 
    ch_filtered = STANDARDISE_AND_FILTER_VARIANTS.out.mutations


}


// if (params.combimets) {
//     // for combimets, additionally assign mutations to clones. 
//     // filter to only mutations in clones on path from trunk to seeds (metastatic trajectory)
//     ch_trees_dir = Channel.fromPath("${params.basedir}/data/phylogeny_angel/trees", type: 'dir')
//     ch_ccfs_dir = Channel.fromPath("${params.basedir}/data/phylogeny_angel/dpclust", type: 'dir')
//     mettraj_clones = file("${params.basedir}/data/phylogeny_angel/met_trajectory_clones.tsv")
//     ch_mutations = ASSIGN_CLONES_COMBIMETS(MERGE_VARIANTS.out.merged, ch_trees_dir, ch_ccfs_dir, mettraj_clones, ch_samplesheet)
// } else {
//     // general
//     ch_mutations = ASSIGN_CLONES_PLACEHOLDER(MERGE_VARIANTS.out.merged)
// }


    // // assign mutations to clones. 
    // // any donor with 2+ samples is expected to have a phylogenetic tree. 
    // ch_trees_dir = Channel.fromPath("${params.basedir}/data/phylogeny_angel/trees", type: 'dir')
    // ch_ccfs_dir = Channel.fromPath("${params.basedir}/data/phylogeny_angel/dpclust", type: 'dir')
    // mettraj_clones = file("${params.basedir}/data/phylogeny_angel/met_trajectory_clones.tsv")
    // ch_assigned = ASSIGN_CLONES(ch_filtered, ch_trees_dir, ch_ccfs_dir, mettraj_clones, ch_samplesheet) 

    // // --- COHORT CONTRASTS --- //
    // // define contrasts 
    // ch_cohorts      = ch_samplesheet | splitCsv(header: true, sep: '\t') | map { row -> row.cohort } | unique
    // ch_contrasts    = ch_cohorts    | combine(ch_data) | combine(ch_samplesheet) | SPLIT_MUTATIONS
    
    // // for each contrast do...
    // ch_contrasts | SINGLE_GENE_ANALYSIS
    // ch_selected     = ch_contrasts | SELECT_GENESETS
    // ch_dif_mutate   = ch_contrasts | join(ch_selected) | DIFFERENTIAL_MUTATION
    // ch_bkg_enrich   = ch_contrasts | join(ch_selected) | BACKGROUND_ENRICHMENT
    // ch_results      = ch_contrasts | join(ch_selected) | join(ch_dif_mutate) | join(ch_bkg_enrich) | SUMMARISE_RESULTS