
include { LOG_SETTINGS } from './modules/misc'

include { EXTRACT_STRUCTVARS } from './modules/extract'
include { EXTRACT_SEQVARS as EXTRACT_SNVS } from './modules/extract'
include { EXTRACT_SEQVARS as EXTRACT_INDELS } from './modules/extract'
include { EXTRACT_CNA } from './modules/extract'
include { MERGE_DONORS } from './modules/extract'

include { HARMONISE_DATA } from './modules/contrast'
include { SPLIT_MUTATIONS } from './modules/contrast'
include { SELECT_GENESETS } from './modules/contrast'
include { DIFFERENTIAL_MUTATION } from './modules/contrast'
include { BACKGROUND_ENRICHMENT } from './modules/contrast'
include { SUMMARISE_RESULTS } from './modules/contrast'

ch_samplesheet  = Channel.fromPath(params.samplesheet) 
genesets        = file(params.genesets)
sizes           = file(params.sizes)
gff             = file(params.refseq_gff_path)
sv_files        = Channel.fromPath("${params.basedir}/data/sv/*.vcf")
snv_files       = Channel.fromPath("${params.basedir}/data/snv/*.vcf")
indel_files     = Channel.fromPath("${params.basedir}/data/indel/*.vcf")
scna_files      = Channel.fromPath("${params.basedir}/data/scna/*.txt")

Channel.fromPath("${params.basedir}/samplesheet.tsv")  
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

all_svs     = sv_files      | map { filepath -> tuple(filepath.simpleName.toString()[0..8], filepath) }
all_snvs    = snv_files     | map { filepath -> tuple(filepath.simpleName.toString()[0..8], filepath) }
all_indels  = indel_files   | map { filepath -> tuple(filepath.simpleName.toString()[0..8], filepath) }
all_scna    = scna_files    | map { filepath -> tuple(filepath.simpleName.toString()[0..8], filepath) }

ch_indels_raw   = samples_ch | join(all_indels) | join(all_scna)
ch_snvs_raw     = samples_ch | join(all_snvs) | join(all_scna)
ch_svs_raw      = samples_ch | join(all_svs)
ch_cna_raw      = samples_ch | join(all_scna)


workflow {

    LOG_SETTINGS()

    // EXTRACT VARIANTS
    ch_svs      = EXTRACT_STRUCTVARS(ch_svs_raw)
    ch_indels   = EXTRACT_INDELS(ch_indels_raw, 'INDEL')
    ch_snvs     = EXTRACT_SNVS(ch_snvs_raw, 'SNV')
    ch_cna      = EXTRACT_CNA(ch_cna_raw, gff)
    MERGE_DONORS(
        ch_svs.collect(),
        ch_cna.collect(),
        ch_snvs.collect(),
        ch_indels.collect()
    )
    ch_mutations = MERGE_DONORS.out.merged
    ch_data = HARMONISE_DATA(ch_mutations, genesets, sizes) 
    
    // TODO merge samples for each donor combimets 
    // TODO filter at later date. 
    // does the CADD feature annotation mirror the VCF ANNOVAR feature annotation? No.
    // ANNOVAR and CADD use different transcripts. 
    // CADD scores should be used with caution. 
    // Eg. PPCG0086a - chr1:9664540 G>T (TMEM201). 
    //   - ANNOVAR says 3_prime_UTR_variant, CADD says intron_variant

    // COHORT CONTRASTS
    // define contrasts 
    ch_cohorts      = ch_samplesheet | splitCsv(header: true, sep: '\t') | map { row -> row.cohort } | unique
    ch_contrasts    = ch_cohorts    | combine(ch_data) | combine(ch_samplesheet) | SPLIT_MUTATIONS
    
    // for each contrast do...
    ch_selected     = ch_contrasts | SELECT_GENESETS
    ch_dif_mutate   = ch_contrasts | join(ch_selected) | DIFFERENTIAL_MUTATION
    ch_bkg_enrich   = ch_contrasts | join(ch_selected) | BACKGROUND_ENRICHMENT
    ch_results      = ch_contrasts | join(ch_selected) | join(ch_dif_mutate) | join(ch_bkg_enrich) | SUMMARISE_RESULTS

}

