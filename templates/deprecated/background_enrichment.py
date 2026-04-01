
import argparse
import pandas as pd 
import numpy as np
from scipy.stats import binomtest
# from statsmodels.stats.multitest import multipletests

def main() -> None:
    args = load_cmdline_args()

    # load data 
    gframe_full, gframe = load_genesets(args.genesets, args.selection)
    sizes = pd.read_csv(args.sizes, sep='\t', header=0)
    muts = pd.read_csv(args.mutations, sep='\t', header=0)
    muts = muts[muts['vtype']!='CNApeak'].copy()
    genelist = generate_genelist(muts, sizes, gframe_full)
    sizes = sizes.set_index('gene')

    # enrichment
    counts_obs = calc_observed(muts, genelist, gframe)
    counts_exp = calc_background(muts, genelist, gframe, sizes)
    counts = counts_obs.copy()
    counts['background'] = counts_exp['background']
    print(counts.head())
        
    # final formatting
    results = counts.copy()
    n_samples = muts['donor'].nunique()
    results['pval'] = results.apply(calc_pval, field='total', n_samples=n_samples, axis=1)
    results['factor'] = results['total'] / results['background']
    results['factor_lvl'] = results['factor'].apply(lambda x: round(x, 1))
    results = results.sort_values(['factor_lvl', 'pval'], ascending=[False, True])
    results = results.drop('factor_lvl', axis=1)
    results = results.reset_index()
    results = results.rename(columns={'index': 'geneset'})
    results = results[[x for x in results.columns if x!='geneset']+['geneset']].copy()
    results.to_csv(args.outfile, sep='\t', index=False, float_format='%.3f')
    print(results.head(20))


###############
### FILE IO ###
###############

def load_genesets(genesets_path: str, selection_path: str):
    # load data 
    gframe = pd.read_csv(genesets_path, sep='\t', header=0)
    sel = pd.read_csv(selection_path, sep='\t', header=0)
    
    # remove Mitochondrial genes
    gframe = gframe[~gframe['gene'].str.startswith('MT-')].copy()
    
    # make new genesets frame by subsetting selected
    all_genesets = sorted(list(gframe['geneset'].unique()))
    sel_genesets = sorted(list(sel[sel['valid']==True]['geneset'].unique()))
    assert len(sel_genesets) != 0
    assert len(set(sel_genesets) - set(all_genesets)) == 0
    sframe = gframe[gframe['geneset'].isin(sel_genesets)].copy()
    return gframe, sframe

def load_sizes(filepath: str) -> pd.DataFrame:
    # load gene length information
    df = pd.read_csv(filepath, sep='\t', header=0)
    df = df.set_index('gene')
    return df

def generate_genelist(muts: pd.DataFrame, sizes: pd.DataFrame, gframe: pd.DataFrame) -> list[str]:
    # load gene length information
    gset_genes = set(gframe['gene'].unique())
    size_genes = set(sizes['gene'].unique())
    mut_genes = set(muts['gene'].unique())
    all_genes = gset_genes | size_genes | mut_genes
    return sorted(list(all_genes))


##################
### ENRICHMENT ###
##################

def gen_dtable(muts: pd.DataFrame, genelist: list[str]) -> pd.DataFrame:
    donors = sorted(list(muts['donor'].unique()))
    dtable = pd.DataFrame(index=genelist)
    donors2genes = muts.groupby('donor')['gene'].agg(set)
    i = 0
    for donor in donors:
        genes = sorted(list(donors2genes[donor]))
        dtable[donor] = 0
        dtable.loc[genes, donor] = 1
        i += 1
        if i % 50 == 0:
            dtable = dtable.copy()
    return dtable

def calc_observed(table: pd.DataFrame, genelist: list[str], gframe: pd.DataFrame) -> pd.DataFrame:
    runs = [
        (['SNV'], 'SNV'), 
        (['INDEL'], 'INDEL'), 
        (['SV'], 'SV'), 
        (['CNA'], 'CNA'), 
        (['SNV', 'INDEL', 'SV', 'CNA'], 'total'), 
    ]
    genesets = sorted(list(gframe['geneset'].unique()))
    gset2genes = gframe.groupby('geneset')['gene'].agg(set)
    counts = pd.DataFrame(index=genesets)

    # calc expected background individually for seqvars/structvars/CNA. 
    for vclasses, mtype in runs:
        df = table[table['vclass'].isin(vclasses)].copy()
        dtable = gen_dtable(df, genelist)
        for gset in genesets:
            genes = sorted(list(gset2genes[gset]))
            dslice = dtable.loc[genes]
            n_donors = dslice.sum().clip(0, 1).sum()
            counts.loc[gset, mtype] = n_donors 
        counts[mtype] = counts[mtype].astype(int)
    assert counts.isna().sum().sum() == 0
    return counts 

def calc_background(table: pd.DataFrame, genelist: list[str], gframe: pd.DataFrame, sizes: pd.DataFrame) -> pd.DataFrame:
    runs = [
        (['SNV', 'INDEL'], 'snv/indel'), 
        (['SV', 'CNA'], 'sv/cna'), 
    ]

    # calc expected background individually for seqvars/structvars/CNA.
    n_samples = table['donor'].nunique()
    expect = pd.DataFrame(index=genelist)
    for vclasses, mtype in runs:
        df = table[table['vclass'].isin(vclasses)].copy()
        dtable = gen_dtable(df, genelist)
        
        # prepare final inputs 
        mutation_matrix = dtable.to_numpy()
        size_field = 'cum_exon_len' if mtype == 'snv/indel' else 'span'
        gene_names = dtable.index.to_list()
        gene_lengths = [int(sizes.loc[gene, size_field]) for gene in genelist]

        # do the calculation
        results = _do_background(mutation_matrix, np.array(gene_lengths), np.array(gene_names))
        results = results.set_index('name')
        expect[f"obs. {mtype}"] = results['observed_mutations']
        expect[f"exp. {mtype}"] = results['expected_mutations']

    exp_cols = [x for x in expect.columns if x.startswith('exp. ')]
    expect['exp. total'] = expect.apply(calc_exp_total, n_samples=n_samples, exp_cols=exp_cols, axis=1)
    # print()
    # print(expect.head(10))

    genesets = sorted(list(gframe['geneset'].unique()))
    gset2genes = gframe.groupby('geneset')['gene'].agg(set)
    counts = pd.DataFrame(index=genesets)
    for gset in genesets:
        genes = sorted(list(gset2genes[gset]))
        counts.loc[gset, 'background'] = expect.loc[genes, 'exp. total'].sum()

    assert counts.isna().sum().sum() == 0
    return counts 

def _do_background(mutation_matrix, gene_lengths, gene_names) -> pd.DataFrame:
    n_genes, n_samples = mutation_matrix.shape
        
    # Per-sample mutation counts
    mutations_per_sample = mutation_matrix.sum(axis=0)
    total_length = np.sum(gene_lengths)

    results = []
    for i in range(n_genes):
        obs = mutation_matrix[i, :].sum()
        
        # Calculate expected mutations accounting for per-sample burden
        expected = 0
        for j in range(n_samples):
            # Given that sample j has mutations_per_sample[j] mutations,
            # what's the expected number in gene i?
            # This follows a hypergeometric-like distribution, but we approximate
            # with: P(gene i hit) ≈ 1 - (1 - gene_length_i/total_length)^mutations_in_sample_j
            
            # Simpler approximation for rare mutations:
            prob_gene_i_in_sample_j = (gene_lengths[i] / total_length) * mutations_per_sample[j]
            expected += min(prob_gene_i_in_sample_j, 1.0)  # Cap at 1 (can't exceed binary)
        
        # Background rate per sample
        gene_background_rate = expected / n_samples
        results.append({
                'gene_idx': i,
                'total_length': gene_lengths[i],
                'observed_mutations': int(obs),
                'expected_mutations': expected,
                'gene_background_rate': gene_background_rate,
            })

    df = pd.DataFrame(results)
    # names
    df['name'] = gene_names
    cols = ['name'] + [col for col in df.columns if col != 'name']
    df = df[cols]
    df = df.drop('gene_idx', axis=1)
    return df 

def calc_exp_total(row: pd.Series, n_samples: int, exp_cols: list[str]) -> float:
    a = row[exp_cols[0]] / n_samples
    for i in range(1, len(exp_cols)):
        b = row[exp_cols[i]] / n_samples
        t = a*b 
        u = b-t 
        a = a+u 
    return a * n_samples

def calc_pval(row: pd.Series, field: str, n_samples: int) -> float:
    obs = int(row[field])
    expected = float(row['background'])
    
    # Background rate per sample
    background_rate = expected / n_samples
    
    # Test observed vs expected
    p_val = binomtest(obs, n_samples, background_rate, alternative='two-sided').pvalue
    return float(p_val)


####################
### CMDLINE ARGS ###
####################


def load_cmdline_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Performs differential mutation analysis between two cohorts.')
    parser.add_argument('--mutations', type=str, required=True, help='Path to mutations (tsv).')
    parser.add_argument('--genesets', type=str, required=True, help='Path to genesets (tsv).')
    parser.add_argument('--selection', type=str, required=True, help='Path to genesets selection (tsv).')
    parser.add_argument('--sizes', type=str, required=True, help='Path to gene sizes (tsv).')
    parser.add_argument('--outfile', type=str, required=True, help='Output file containing results (tsv).')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    main()