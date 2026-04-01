
GMT_ONCOGENIC='/home/grace/work/PPCG_DifferentialGenesetMutation/data/genesets/c6.all.v2026.1.Hs.symbols.gmt'
GMT_GAVISH='/home/grace/work/PPCG_DifferentialGenesetMutation/data/genesets/c4.3ca.v2026.1.Hs.symbols.gmt'
GMT_KEGG='/home/grace/work/PPCG_DifferentialGenesetMutation/data/genesets/c2.cp.kegg_medicus.v2026.1.Hs.symbols.gmt'
GMT_GOBP='/home/grace/work/PPCG_DifferentialGenesetMutation/data/genesets/c5.go.bp.v2026.1.Hs.symbols.gmt'
GMT_REACTOME='/home/grace/work/PPCG_DifferentialGenesetMutation/data/genesets/c2.cp.reactome.v2026.1.Hs.symbols.gmt'

MUTATIONS='/home/grace/work/PPCG_DifferentialGenesetMutation/data/manual/mutations.assigned.010426.tsv'
SIZES='/home/grace/work/PPCG_DifferentialGenesetMutation/outputs/alldonors_01042026/variant_processing/sizes.tsv'
SHEET='/home/grace/work/PPCG_DifferentialGenesetMutation/samplesheet.angel.alldonors.tsv'
CLONE_META='/home/grace/work/PPCG_DifferentialGenesetMutation/data/angel/26.03.2026/metastatic_clones_with_ancestors.tsv'
# outdir='/home/grace/work/PPCG_DifferentialGenesetMutation/results/latency'
# python run_latency.py --run-id Top50Genes --outdir $outdir --genes
# python run_latency.py --run-id Top20Genes --outdir $outdir --genes
# python run_latency.py --run-id KEGG --outdir $outdir --gmt $GMT_KEGG
# python run_latency.py --run-id GAVISH --outdir $outdir --gmt $GMT_GAVISH
# python run_latency.py --run-id ONCOGENIC --outdir $outdir --gmt $GMT_ONCOGENIC
# python run_latency.py --run-id REACTOME --outdir $outdir --gmt $GMT_REACTOME
# python run_latency.py --run-id GOBP --outdir $outdir --gmt $GMT_GOBP


outdir='/home/grace/work/PPCG_DifferentialGenesetMutation/test/metastasis'
# python run_metastasis.py --run-id KEGG --gmt $GMT_KEGG --muts $MUTATIONS --sizes $SIZES --sheet $SHEET --clone-meta $CLONE_META --outdir $outdir 
python run_metastasis.py --run-id REACTOME --gmt $GMT_REACTOME --muts $MUTATIONS --sizes $SIZES --sheet $SHEET --clone-meta $CLONE_META --outdir $outdir 

