
GMT_ONCOGENIC='/home/grace/work/PPCG_DifferentialGenesetMutation/data/genesets/c6.all.v2026.1.Hs.symbols.gmt'
GMT_GAVISH='/home/grace/work/PPCG_DifferentialGenesetMutation/data/genesets/c4.3ca.v2026.1.Hs.symbols.gmt'
GMT_KEGG='/home/grace/work/PPCG_DifferentialGenesetMutation/data/genesets/c2.cp.kegg_medicus.v2026.1.Hs.symbols.gmt'
GMT_GOBP='/home/grace/work/PPCG_DifferentialGenesetMutation/data/genesets/c5.go.bp.v2026.1.Hs.symbols.gmt'
GMT_REACTOME='/home/grace/work/PPCG_DifferentialGenesetMutation/data/genesets/c2.cp.reactome.v2026.1.Hs.symbols.gmt'

outdir='/home/grace/work/PPCG_DifferentialGenesetMutation/results/latency'
# python run_latency.py --run-id Top50Genes --outdir $outdir --genes
python run_latency.py --run-id Top20Genes --outdir $outdir --genes
python run_latency.py --run-id KEGG --outdir $outdir --gmt $GMT_KEGG
# python run_latency.py --run-id GAVISH --outdir $outdir --gmt $GMT_GAVISH
# python run_latency.py --run-id ONCOGENIC --outdir $outdir --gmt $GMT_ONCOGENIC
# python run_latency.py --run-id REACTOME --outdir $outdir --gmt $GMT_REACTOME
# python run_latency.py --run-id GOBP --outdir $outdir --gmt $GMT_GOBP

