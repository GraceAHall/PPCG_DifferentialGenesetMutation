

import argparse 
import numpy as np
import pandas as pd 
import statsmodels.api as sm

from scipy import stats
from numpy.linalg import LinAlgError

import warnings 
warnings.filterwarnings('ignore')


def main() -> None:
    args = load_cmdline_args()
    table, gframe = load_data(args)
    matrix = gen_matrix(table, gframe)
    rframe = evaluate(matrix, gframe)
    rframe.to_csv(args.outfile, sep='\t', index=False, float_format='%.3f')

def load_data(args: argparse.Namespace):
    # mutations
    pmuts = pd.read_csv(args.posmuts, sep='\t', header=0)
    nmuts = pd.read_csv(args.negmuts, sep='\t', header=0)
    pmuts['cls'] = 'positive'
    nmuts['cls'] = 'negative'
    table = pd.concat([pmuts, nmuts], ignore_index=True)
    table = table[table['vtype']!='CNApeak'].copy()

    # genesets, filtering
    gframe = pd.read_csv(args.genesets, sep='\t', header=0)
    sframe = pd.read_csv(args.selection, sep='\t', header=0)
    sel_genesets = sorted(list(sframe[sframe['valid']==True]['geneset'].unique()))
    all_genesets = sorted(list(gframe['geneset'].unique()))
    assert len(sel_genesets) != 0
    assert len(set(sel_genesets) - set(all_genesets)) == 0
    gframe = gframe[gframe['geneset'].isin(sel_genesets)].copy()
    return table, gframe

def gen_matrix(table: pd.DataFrame, gframe: pd.DataFrame) -> pd.DataFrame:
    # init dataframe
    all_donors = sorted(list(table['donor'].unique()))
    df = pd.DataFrame(index=all_donors)

    # annotate membership & TMB
    df['cls'] = table.drop_duplicates('donor').set_index('donor')['cls']
    df['TMB'] = table.groupby('donor')['gene'].nunique() #? coords or genes....hmmm
    df = df.reset_index()
    df = df.rename(columns={'index': 'donor'})
    df = df[['cls', 'donor', 'TMB']].copy()
    df = df.sort_values('cls', ascending=False)

    # for each geneset, add column annotating whether the donor has a mutation or not. 
    gene2donors = table.groupby('gene')['donor'].agg(set).to_dict()
    gset2genes = gframe.groupby('geneset')['gene'].agg(set).to_dict()
    i = 0
    for gset, genes in gset2genes.items():
        g2d_sub = {k: v for k, v in gene2donors.items() if k in genes}
        donors = set.union(*list(g2d_sub.values()))
        df[gset] = df['donor'].apply(lambda x: 1 if x in donors else 0)
        i += 1
        if i % 10 == 0:
            df = df.copy()

    return df

def evaluate(df: pd.DataFrame, gframe: pd.DataFrame) -> pd.DataFrame:
    all_genesets = sorted(list(gframe['geneset'].unique()))
    data = []
    for gset in all_genesets:
        counts = df.groupby('cls')[gset].sum()
        fish_pval = _calc_fisher(df, gset)
        logit_pval, converged = _calc_logit(df, gset)
        data.append((gset, counts['positive'], counts['negative'], fish_pval, logit_pval, converged))

    rframe = pd.DataFrame.from_records(data, columns=['geneset', 'posclass donors', 'negclass donors', 'fisher_pval', 'logit_pval', 'logit_converged'])
    rframe = rframe.sort_values('logit_pval')
    return rframe

def _calc_fisher(df: pd.DataFrame, group: str) -> float:
    a = ((df[group]==1) & (df['cls']=='positive')).sum()
    b = ((df[group]==1) & (df['cls']=='negative')).sum()
    c = ((df[group]==0) & (df['cls']=='positive')).sum()
    d = ((df[group]==0) & (df['cls']=='negative')).sum()
    res = stats.fisher_exact([[a, b], [c, d]], alternative='two-sided')
    return float(res.pvalue)  # type: ignore

def _calc_logit(df: pd.DataFrame, group: str):
    y = df[group]
    X = pd.DataFrame({
        "cls": (df["cls"] == "positive").astype(int),
        "log_TMB": np.log1p(df["TMB"]) # log(1 + TMB)
    })
    X = sm.add_constant(X)
    try:
        model = sm.Logit(y, X).fit(disp=0)
        return float(model.pvalues['cls']), model.mle_retvals['converged']
    except LinAlgError:
        return 1.0, False


def load_cmdline_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Performs differential mutation analysis between two cohorts.')
    parser.add_argument('--posmuts', type=str, required=True, help='Path to positive class mutations (tsv).')
    parser.add_argument('--negmuts', type=str, required=True, help='Path to negative class mutations (tsv).')
    parser.add_argument('--genesets', type=str, required=True, help='Path to genesets file (after harmonisation).')
    parser.add_argument('--selection', type=str, required=True, help='Path to genesets selection tsv file.')
    parser.add_argument('--outfile', type=str, required=True, help='Output file containing results (tsv).')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    main()