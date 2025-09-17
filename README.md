
# PPCG Gene Set Mutation Analysis Pipeline

## Setup

**Data and software** 

This pipeline has been tested on nextflow version 24.10.6 and Python 3.11. 

Multiple pieces of data are required. These should be in a `data/` dir within the pipeline directory. 

Before you start, ensure your workspace has the following structure. 

<pre>
├── README.md
├── data
│   ├── genesets    # gene set lists (only one is used)
│   ├── indel       # Raw PPCG INDEL VCFs
│   ├── meta.csv    
│   ├── scna        # Battenberg SCNA segments
│   ├── snv         # Raw PPCG SNV VCFs
│   ├── supporting  # Supporting files
│   └── sv          # SnpEff annotated SV VCFs
├── main.nf
├── modules
├── nextflow.config
├── requirements.txt
├── samplesheet.tsv
└── templates       # python scripts used in pipeline
</pre>


**Python Environment**


Before running the pipeline, please follow these steps. 

1. Ensure python and nextflow have equal or later version. 
- `python -V` should be 3.11 or above. 
- `nextflow -v` should be 24.10.6 or above. 

2. Create a virtual environment using python `venv`<br>
`python -m venv venv`

3. Activate the virtual environment<br>
`source venv/bin/activate`

4. Install python requirements using pip<br>
`pip install -r requirements.txt`


## Pipeline Summary

This pipeline performs gene set mutation enrichment between two cohorts. All settings for each process are available in the `nextflow.config` file and should be tailored to your analysis. 

Main steps:
1. Log settings 
2. Extract variants 
    - Extract SNV/INDEL/SV/CNA variants from each donor
    - Merge all variants into single file & provide summary.
    - Harmonise gene symbols
3. Gene set mutation enrichment
    - Select gene sets for analysis
    - Differential mutation between cohorts 
    - Background enrichment within each cohort
    - Summarise results

**Logging settings**

For each workflow execution, please provide a meaningful `params.run_id` name. 
This should reflect the settings you have altered in the `nextflow.config` file. 
All settings are logged in `outputs/<run_id>/pipeline_info/settings.txt` for later reference. 
This is important so you can keep track of different outputs when settings are changed. 

**SNV extraction**

SNV and INDEL variants are extracted from each donor's raw VCF file. ANNOVAR annotations are used to define the impacted feature and function. 

By default, all coding region variants are extracted (except synonymous variants). UTR, splice, noncoding variants, and upstream/downstream variants can be optionally enabled. Intron variants are not extracted. 

The following `params.seqvar` settings can be customised:
- `allow_utr`: Include 3' and 5' UTR variants.
- `allow_splice`: Include splicing variants. 
- `allow_noncoding`: Include noncoding variants. 
- `allow_intergenic_dist`: Include upsteam/downstream variants within \<dist> of gene.

**SV extraction**

SVs are extracted from each donor's SnpEff VCF file. SnpEff was run on the raw VCF file to provide feature impact annotations. Any variant with a 'HIGH' SnpEff impact annotation is extracted. These are fusions and transcript ablations. No settings are currently available for customisation. 

**CNA extraction**

Copy number alterations are extracted from Battenberg subclonal copy number segments. The type of alterations (AMP, shallow DEL, deep DEL) to extract are available as options. Thresholds to call each variant type depends on LOH and WGD status. The exact logic can be see in the three `CNAEventValidator` classes within `templates/extract_cna.py`.

The following `params.cnavar` settings can be customised:
- `min_span`: Minimum affected proportion of the gene to call a variant. 
- `allow_amp`: Include amplifications.
- `allow_deep_del`: Include deep deletions (total copy number == 0)
- `allow_shallow_del`: Include shallow deletions (total copy number <= 1 if WGD, else == 0).
- `allow_subclonal`: Include Battenberg segments where two states are predicted (subclonal).
- `allow_noncoding`: Include non-coding genes. 
- `allow_chrX`: Include chromosome X variants. Chromosome Y not available. 

**Variant merging**

After all variants for all donors have been extracted, individual files are merged into a single `merged.tsv` file representing all donors.  This one file contains all variant types for all donors in `samplesheet.tsv`. 

Summary statistics about extracted variants are available in `merged.log`. This file documents:
- Top mutated genes
- Number of donors mutated per broad class (SNV, INDEL, SV, CNA)
- Number of donors mutated per subclasses <br>eg `missense_variant`, `splice_site_variant` for SNVs,  `bidirectional_gene_fusion` for SVs etc
- The above 2 points, but stratified by gene. <br>eg PTEN
    - 68 donors have CNA
    - 13 donors have INDEL
    - 13 donors have SNV
    - 99 donors have SV

**Harmonising gene symbols**

Gene symbols are standardised across sources before downstream analysis. This ***aims*** to ensure no genes are slipping through the cracks, causing us to miss gene sets which would otherwise be significant. 

Symbols used in this pipeline are derived from three sources. 
- PPCG mutations: Genes annotated by ANNOVAR / SnpEff.
- RefSeq: Genes used in RefSeq GRCh37 GFF file.
- GSEA: Genes used by MSigDB gene sets (eg C5 GO:BP).

These sometimes conflict as symbols can have aliases or be updated / deprecated over time. For example if a given gene has a different symbol in the MSigDB gene set vs the extracted mutation data (ANNOVAR/SnpEff) we may  accidently ignore the gene. Ideally, the pipeline should use identifiers rather than symbols (will implement eventually). For the time being, these need to be harmonised across all sources.  

**Selecting gene sets**

Gene sets are selected before downstream analysis. This is to reduce the number of gene sets being assessed. The choice of parameter values has a huge impact on results. 

There are two provided MSigDB gene set lists in the `data/genesets` folder. 
- `c2.cp.kegg_medicus.v2024.1.Hs.symbols.tsv` (658 gene sets) <br>Canonical Pathways gene sets derived from the KEGG MEDICUS pathway database.
- `c5.go.bp.v2024.1.Hs.symbols.tsv` (7583 gene sets) <br>Gene sets derived from the GO Biological Process ontology.

The following `params.select_genesets` settings can be customised:
- `min_genes`: Minimum number of genes in gene set. 
- `max_genes`: Maximum number of genes in gene set.
- `min_genes_hitprop`: Minimum proportion of ***genes*** with mutation in gene set.<br>eg if there are 10 genes in the gene set, a value of `0.2` means 2 separate genes must be mutated in the cohort. 
- `min_cohort_hitprop`: Minimum proportion of ***donors*** in cohort with mutation in gene set.<br>eg if there are 100 donors, a value of `0.2` means 20 donors must have a mutation in the gene set. 
- `max_gene_reliance`: The maximum reliance on a single gene in the gene set.<br>eg. `0.5` would remove gene sets where a single gene accounts for 50% of all donor hits. 

The `max_gene_reliance` parameter is especially useful as can be used to remove gene sets dominated by mutations in a single gene, eg for PTEN, TP53, ERG etc. 

**Differential mutation**

Differential mutation is assessed using logistic regression. This performs the same role as a Fisher test, but includes the mutational burden of each donor as a covariate. I.e. donors with high mutational burden are weighted lower than those with low burden.  

Two groups are required - a 'positive' group and a 'negative' group. Eg. for CIN.Cluster-2 in the mutational processes group, the *positive* group is donors marked as CIN.Cluster-2, and the *negative* class is donors in other clusters. 

**Background enrichment**

Enrichment above background is a secondary test to ensure that mutations in a given gene set are observed more frequently than pure chance. Bootstrapping is used to permute gene mutation assignments per donor (preserves mutational burden), then a p-value is calculated as the number of bootstrap samples where the number of affected donors is equal or greater than the 'actual' observed data. 

Permutation method during bootstrap sample generation depends on variant class:
- SNVs/INDELs: Selection probability for a given gene is weighted by cumulative exon length. 
- SVs: Selection probability for a given gene is weighted by total gene span. 
- CNA: No probabilities, just randomly select a gene. 

The background enrichment process has the longest runtime by far and scales with the value for `budget`. The `params.background_enrichment.budget` setting controls this value by specifying how many bootstrap samples are performed. A higher budget is a more accurate estimate of the background mutational rate, at the cost of runtime. A minimum of 100 bootstrap samples (budget=100) is recommended for robustness, but feel free to reduce this (eg budget=10) if doing test runs to save time. 

**Summarising results**

The final process collates results into a single table and generates OncoPrints. 

The p-values are unadjusted. The resposibility for multiple testing adjustment falls on the user. This choice has been made as different users may have different methods for adjustment, and may use different schemes when considering the number of 'tests' which have been performed. For example, a reasonable process may be to first remove gene sets not enriched above background within the cohort before considering the differential mutation (logistic regression) results, thus reducing the amount of gene sets being 'tested'. 

OncoPrints are generated using PyComplexHeatmap so that individual mutations can be viewed per gene set. Only the top 'significant' gene sets have an oncoprint generated. The number of OncoPrints to create can be changed using the `params.summarise_results.max_plots` parameter in `nextflow.config`. For example a value of `10` will create oncoprints for the top 10 most significant genesets. 

The final results table has the following fields:
- *index*: Name of the gene set.
- *obs_donors*: Number of donors with a mutation in the gene set.
- *obs_proportion*: Proportion of cohort with a mutation in the gene set.
- *bkg_donors*: Estimated background number of donors expected to have a mutation by pure chance. 
- *bkg_proportion*: Estimated background proportion of cohort expected to have a mutation by pure chance. 
- *obs_logit_pval*: The p-value for trained logistic regression model. 
- *obs_factor*: Ratio of positive vs negative cohort mutated (proportion of each cohort).
- *bkg_enrich_factor*: Ratio of positive vs negative cohort mutation above background (proportion of each cohort). 

The `obs_donors`, `obs_proportion`, `bkg_donors` and `bkg_proportion` fields are represented as `<positive> / <negative>`, where *positive* refers to the treatment group, and *negative* refers to the control group. This is to avoid having too many columns in the tsv file. Eg. for CIN.Cluster-2 in the mutational processes group, the *positive* group is donors marked as CIN.Cluster-2, and the *negative* class is donors in other clusters. 


