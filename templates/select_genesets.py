
import argparse 
import pandas as pd 
from typing import Tuple
from collections import defaultdict


def main() -> None:
    args = load_cmdline_args()

    pmuts = pd.read_csv(args.posmuts, sep='\t', header=0)
    nmuts = pd.read_csv(args.negmuts, sep='\t', header=0)
    gframe = pd.read_csv(args.genesets, sep='\t', header=0)
    pmuts = pmuts[pmuts['vtype']!='CNApeak'].copy()
    nmuts = nmuts[nmuts['vtype']!='CNApeak'].copy()

    sizes               = get_size(gframe)
    n_donors_hit        = get_donors_hit(gframe, pmuts)
    n_genes_hit         = get_genes_hit(gframe, pmuts)
    topgenes, reliances = get_topgene_reliance(gframe, pmuts)

    sframe = pd.DataFrame(index=sorted(list(gframe['geneset'].unique())))
    sframe['size'] = sizes
    sframe['donors_hit'] = n_donors_hit
    sframe['genes_hit'] = n_genes_hit
    sframe['topgene'] = topgenes
    sframe['reliance'] = reliances

    n_donors = pmuts['donor'].nunique()
    sframe['n_genes_ok'] = (sframe['size']>=args.min_genes) & (sframe['size']<=args.max_genes)
    sframe['gene_hitprop_ok'] = sframe['genes_hit'] / sframe['size'] >= args.min_genes_hitprop
    sframe['cohort_hitprop_ok'] = sframe['donors_hit'] / n_donors >= args.min_cohort_hitprop
    sframe['gene_reliance_ok'] = sframe['reliance'] <= args.max_gene_reliance
    sframe['valid'] = sframe[['n_genes_ok', 'gene_hitprop_ok', 'cohort_hitprop_ok', 'gene_reliance_ok']].all(axis=1)

    sframe = sframe.reset_index()
    sframe = sframe.rename(columns={'index': 'geneset'})
    sframe = sframe[['geneset', 'valid'] + [x for x in sframe.columns if x not in ['valid', 'geneset']]]
    sframe = sframe.sort_values(['valid', 'donors_hit'], ascending=[False, False])
    sframe.to_csv(args.outfile, sep='\t', index=False, float_format='%.2f')


def get_size(df: pd.DataFrame) -> pd.Series:
    return df.groupby('geneset')['gene'].nunique()

def get_donors_hit(df: pd.DataFrame, muts: pd.DataFrame) -> pd.Series:
    gene2gset = df.groupby('gene')['geneset'].agg(set).to_dict()
    gset2donors = defaultdict(set)
    for rec in muts.itertuples():
        if rec.gene in gene2gset:
            for gset in gene2gset[rec.gene]:
                gset2donors[gset].add(rec.donor)
    
    data = {}
    all_genesets = sorted(list(df['geneset'].unique()))
    for gset in all_genesets:
        data[gset] = len(gset2donors[gset])
    return pd.Series(data)

def get_genes_hit(df: pd.DataFrame, muts: pd.DataFrame) -> pd.Series:
    mut_genes = set(muts['gene'].unique())
    gset2genes = df.groupby('geneset')['gene'].agg(set).to_dict()
    data = {}
    for gset, genes in gset2genes.items():
        shared = genes & mut_genes
        data[gset] = len(shared)
    return pd.Series(data)

def get_topgene_reliance(df: pd.DataFrame, muts: pd.DataFrame) -> Tuple[pd.Series, pd.Series]:
    gset2genes = df.groupby('geneset')['gene'].agg(set).to_dict()
    gene2donors = muts.groupby('gene')['donor'].agg(set)
    data_topgene = {}
    data_reliance = {}
    for gset, genes in gset2genes.items():
        g2d_sub = {k: v for k, v in gene2donors.items() if k in genes}
        if len(g2d_sub) < 2:
            data_topgene[gset] = '.' 
            data_topgene[gset] = 1.0 
            continue
        topgene = max(g2d_sub, key=lambda k: len(g2d_sub[k]))
        n_donors = len(set.union(*list(g2d_sub.values())))
        g2d_sub_notop = {k: v for k, v in g2d_sub.items() if k!=topgene}
        n_donors_notop = len(set.union(*list(g2d_sub_notop.values())))
        data_topgene[gset] = topgene
        data_reliance[gset] = 1 - (n_donors_notop/n_donors)
        
    return pd.Series(data_topgene), pd.Series(data_reliance)


def load_cmdline_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Marks genesets as valid or not for downstream analysis.')
    parser.add_argument('--posmuts', type=str, required=True, help='Path to positive class mutations (tsv).')
    parser.add_argument('--negmuts', type=str, required=True, help='Path to negative class mutations (tsv).')
    parser.add_argument('--genesets', type=str, required=True, help='Path to genesets file (after harmonisation).')
    parser.add_argument('--min-genes', type=int, default=10, help='Minimum number of genes in geneset.')
    parser.add_argument('--max-genes', type=int, default=50, help='Maximum number of genes in geneset.')
    parser.add_argument('--min-genes-hitprop', type=float, default=0.1, help='Minimum proportion of hit genes in geneset.')
    parser.add_argument('--min-cohort-hitprop', type=float, default=0.15, help='Minimum proportion of hit donors in positive cohort.')
    parser.add_argument('--max-gene-reliance', type=float, default=0.5, help='The maximum reliance on a single gene in the geneset.\neg 0.5 would remove genesets where a single gene accounts for 50% of all hits.')
    parser.add_argument('--outfile', type=str, required=True, help='Output file containing genesets, stats and verdict (tsv).')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    main()