
import argparse
import pandas as pd 
import numpy as np
from scipy.stats import binomtest
# from statsmodels.stats.multitest import multipletests


def main() -> None:
    args = load_cmdline_args()

    # load genesets
    _, gframe = load_genesets(args.genesets, args.selection)
    all_genesets = sorted(gframe['geneset'].unique())

    # load gene sizes
    sizes = load_sizes(args.sizes, gframe, all_genesets)

    # load mutations
    muts = pd.read_csv(args.mutations, sep='\t', header=0)
    muts = muts[muts['vtype']!='CNApeak'].copy()

    # run enrichment
    results = run_enrichment(muts, all_genesets, gframe, sizes)
    results['factor_lvl'] = results['factor'].apply(lambda x: round(x))
    results = results.sort_values(['pval', 'factor_lvl'], ascending=[True, False])
    results = results.drop('factor_lvl', axis=1)
    results = results.reset_index()
    results = results.rename(columns={'index': 'geneset'})
    results = results[[x for x in results.columns if x!='geneset']+['geneset']].copy()
    results.to_csv(args.outfile, sep='\t', index=False, float_format='%.3f')

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

def load_sizes(filepath: str, gframe: pd.DataFrame, all_genesets: list[str]) -> pd.DataFrame:
    # load gene length information
    df = pd.read_csv(filepath, sep='\t', header=0)
    df = df.set_index('gene')

    # calculate geneset sizes
    # (cumulative length of each geneset for different variant types)
    sizes = pd.DataFrame(index=all_genesets, columns=['seqvar', 'structvar', 'CNA'])
    for geneset in all_genesets:
        genes = sorted(gframe[gframe['geneset']==geneset]['gene'].unique())
        sizes.loc[geneset, 'seqvar'] = int(df.loc[genes]['cum_exon_len'].sum())
        sizes.loc[geneset, 'structvar'] = int(df.loc[genes]['span'].sum())
        sizes.loc[geneset, 'CNA'] = int(df.loc[genes]['span'].sum())
    return sizes

def run_enrichment(table: pd.DataFrame, all_genesets: list[str], gframe: pd.DataFrame, size_frame: pd.DataFrame) -> pd.DataFrame:

    # observed counts
    counts = calc_observed(table, all_genesets, gframe)
    print()
    print(counts.head(20))
    
    # background counts
    runs = [
        (['SNV', 'INDEL'], 'seqvar'), 
        (['SV'], 'structvar'), 
        (['CNA'], 'CNA'), 
    ]

    # calc expected background individually for seqvars/structvars/CNA. 
    all_donors = sorted(table['donor'].unique())
    for vclasses, mtype in runs:

        # generate binary matrix        
        df = table[table['vclass'].isin(vclasses)].copy()
        dtable = pd.DataFrame(data=0, index=all_genesets, columns=all_donors)
        for donor in all_donors:
            genes = sorted(list(df[df['donor']==donor]['gene'].unique()))
            gsets = sorted(list(gframe[gframe['gene'].isin(genes)]['geneset'].unique()))
            dtable.loc[gsets, donor] = 1
        
        # prepare final inputs 
        mutation_matrix = dtable.to_numpy()
        gset_names = dtable.index.to_list()
        gset_lengths = [int(size_frame.loc[gset, mtype]) for gset in gset_names]

        # do the calculation
        results = calc_background(mutation_matrix, np.array(gset_lengths), np.array(gset_names))
        results = results.set_index('name')
        # counts[f"obs_{mtype}"] = results['observed_mutations']
        counts[f"exp_{mtype}"] = results['expected_mutations']

    # combine results from each variant class
    n_samples = len(all_donors)
    print()
    print(n_samples)
    print(counts.sort_values('obs_total', ascending=False).head(10))
    exp_cols = [x for x in counts.columns if x.startswith('exp_')]
    print()
    print(exp_cols)
    counts['exp_total'] = counts.apply(calc_exp_total, n_samples=n_samples, exp_cols=exp_cols, axis=1)
    counts['pval'] = counts.apply(calc_pval, n_samples=n_samples, axis=1)

    # calculate fold enrichment 
    mask = counts['exp_total']>0
    counts.loc[mask, 'factor'] = counts.loc[mask, 'obs_total'] / counts.loc[mask, 'exp_total']
    counts.loc[~mask, 'factor'] = np.inf
    counts = counts.sort_values('pval')

    # return
    return counts

def calc_exp_total(row: pd.Series, n_samples: int, exp_cols: list[str]) -> float:
    a = row[exp_cols[0]] / n_samples
    for i in range(1, len(exp_cols)):
        b = row[exp_cols[i]] / n_samples
        t = a*b 
        u = b-t 
        a = a+u 
    return a * n_samples

def calc_observed(table: pd.DataFrame, all_genesets: list[str], gframe: pd.DataFrame) -> pd.DataFrame:
    runs = [
        (['SNV', 'INDEL'], 'seqvar'), 
        (['SV'], 'structvar'), 
        (['CNA'], 'CNA'), 
        (['SNV', 'INDEL', 'SV', 'CNA'], 'total'), 
    ]
    counts = pd.DataFrame(index=all_genesets)

    # calc expected background individually for seqvars/structvars/CNA. 
    for vclasses, mtype in runs:
        df = table[table['vclass'].isin(vclasses)].copy()
        for geneset in all_genesets:
            genes = list(gframe[gframe['geneset']==geneset]['gene'].unique())
            n_donors = df[df['gene'].isin(genes)]['donor'].nunique()
            counts.loc[geneset, f"obs_{mtype}"] = n_donors 
        counts[f"obs_{mtype}"] = counts[f"obs_{mtype}"].astype(int)
    assert counts.isna().sum().sum() == 0
    return counts 

def calc_background(mutation_matrix, gene_lengths, gene_names) -> pd.DataFrame:
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

def calc_pval(row: pd.Series, n_samples: int) -> float:
    obs = int(row['obs_total'])
    expected = float(row['exp_total'])
    
    # Background rate per sample
    background_rate = expected / n_samples
    
    # Test observed vs expected
    p_val = binomtest(obs, n_samples, background_rate, alternative='two-sided').pvalue
    return float(p_val)



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