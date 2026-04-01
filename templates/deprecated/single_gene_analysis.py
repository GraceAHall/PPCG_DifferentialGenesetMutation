
import sys
import argparse 
import numpy as np
import pandas as pd 
import statsmodels.api as sm

from scipy import stats
from numpy.linalg import LinAlgError

import warnings 
warnings.filterwarnings('ignore')

pd.options.display.float_format = "{:,.3f}".format
pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 200)


def main() -> None:
    args = load_cmdline_args()
    table = load_data(args)
    matrix = gen_matrix(table)
    rframe = evaluate(matrix)
    log_summary(args, table, rframe)
    rframe.to_csv(f"{args.outfile}.tsv", sep='\t', index=False, float_format='%.3f')

def load_data(args: argparse.Namespace) -> pd.DataFrame:
    # mutations
    pmuts = pd.read_csv(args.posmuts, sep='\t', header=0)
    nmuts = pd.read_csv(args.negmuts, sep='\t', header=0)
    pmuts['cls'] = 'positive'
    nmuts['cls'] = 'negative'
    df = pd.concat([pmuts, nmuts], ignore_index=True)
    df = df[df['vtype']!='CNApeak'].copy()
    return df

def gen_matrix(table: pd.DataFrame) -> pd.DataFrame:
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
    i = 0
    for gene, donors in table.groupby('gene')['donor'].agg(set).to_dict().items():
        df[gene] = df['donor'].apply(lambda x: 1 if x in donors else 0)
        i += 1
        if i % 50 == 0:
            df = df.copy()
    return df

def evaluate(df: pd.DataFrame) -> pd.DataFrame:
    all_genes = df.columns.to_list()[3:]
    data = []
    for gene in all_genes:
        counts = df.groupby('cls')[gene].sum()
        fish_pval = _calc_fisher(df, gene)
        logit_pval, converged = _calc_logit(df, gene)
        data.append((gene, counts['positive'], counts['negative'], fish_pval, logit_pval, converged))

    rframe = pd.DataFrame.from_records(data, columns=['gene', 'posclass donors', 'negclass donors', 'fisher_pval', 'logit_pval', 'logit_converged'])
    rframe['fisher_pval'] = rframe['fisher_pval'].fillna(1.0)
    rframe['logit_pval'] = rframe['logit_pval'].fillna(1.0)
    rframe['fisher_pFDR'] = stats.false_discovery_control(rframe['fisher_pval'].to_list())
    rframe['logit_pFDR'] = stats.false_discovery_control(rframe['logit_pval'].to_list())
    rframe = rframe.sort_values('fisher_pval')
    return rframe

def _calc_fisher(df: pd.DataFrame, gene: str) -> float:
    a = ((df[gene]==1) & (df['cls']=='positive')).sum()
    b = ((df[gene]==1) & (df['cls']=='negative')).sum()
    c = ((df[gene]==0) & (df['cls']=='positive')).sum()
    d = ((df[gene]==0) & (df['cls']=='negative')).sum()
    res = stats.fisher_exact([[a, b], [c, d]], alternative='two-sided')
    return float(res.pvalue)  # type: ignore

def _calc_logit(df: pd.DataFrame, gene: str):
    y = df[gene]
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

def log_summary(args: argparse.Namespace, df: pd.DataFrame, rframe: pd.DataFrame) -> None:
    original_stdout = sys.stdout

    # open logfile, write all print statements to file (stdout capture)
    with open(f"{args.outfile}.log", 'w') as f:
        sys.stdout = f

        print('------------------')
        print('--- Basic Info ---')
        print('------------------')
        print()
        print(f"--- Donors ---")
        print(f"Total:    {df['donor'].nunique()}")
        print(f"Positive: {df[df['cls']=='positive']['donor'].nunique()}")
        print(f"Negative: {df[df['cls']=='negative']['donor'].nunique()}")
        print()
        print(f"--- Genes ---")
        print(f"Total:    {df['gene'].nunique()}")
        print(f"Positive: {df[df['cls']=='positive']['gene'].nunique()}")
        print(f"Negative: {df[df['cls']=='negative']['gene'].nunique()}")
        print()
        print(f"--- Variants ---")
        print(f"Total:    {df['coords'].nunique()}")
        print(f"Positive: {df[df['cls']=='positive']['coords'].nunique()}")
        print(f"Negative: {df[df['cls']=='negative']['coords'].nunique()}")

        print()
        print('-----------------------------')
        print('--- Variant Class Summary ---')
        print('-----------------------------')
        for clsmem in ['positive', 'negative']:
            print()
            print(f'[{clsmem.capitalize()}]')
            dfslice = df[df['cls']==clsmem]
            ndonors = dfslice['donor'].nunique()
            counts = pd.DataFrame(index=sorted(list(dfslice['vclass'].unique())))
            counts['variants'] = dfslice['vclass'].value_counts()
            counts['genes'] = dfslice.groupby('vclass')['gene'].nunique()
            counts['donors'] = dfslice.groupby('vclass')['donor'].nunique()
            counts['donors prop.'] = counts['donors'] / ndonors
            print(counts)

        print()
        print('----------------------------------')
        print('--- Variant Annotation Summary ---')
        print('----------------------------------')
        for clsmem in ['positive', 'negative']:
            print()
            print(f'[{clsmem.capitalize()}]')
            dfslice = df[df['cls']==clsmem]
            ndonors = dfslice['donor'].nunique()
            counts = pd.DataFrame(index=sorted(list(dfslice['annotation'].unique())))
            counts['variants'] = dfslice['annotation'].value_counts()
            counts['genes'] = dfslice.groupby('annotation')['gene'].nunique()
            counts['donors'] = dfslice.groupby('annotation')['donor'].nunique()
            counts['donors prop.'] = counts['donors'] / ndonors
            counts = counts.sort_values('variants', ascending=False)
            print(counts)

        print()
        print('-----------------------------------')
        print('--- Gene Summary (top 20 genes) ---')
        print('-----------------------------------')
        for clsmem in ['positive', 'negative']:
            print()
            print(f'[{clsmem.capitalize()}]')
            dfslice = df[df['cls']==clsmem]
            ndonors = dfslice['donor'].nunique()
            counts = pd.DataFrame(index=sorted(list(dfslice['gene'].unique())))
            counts['sv'] = dfslice[dfslice['vclass']=='SV']['gene'].value_counts()
            counts['snvs'] = dfslice[dfslice['vclass']=='SNV']['gene'].value_counts()
            counts['indels'] = dfslice[dfslice['vclass']=='INDEL']['gene'].value_counts()
            counts['cna'] = dfslice[dfslice['vclass']=='CNA']['gene'].value_counts()
            counts['total'] = dfslice['gene'].value_counts()
            counts['donors'] = dfslice.groupby('gene')['donor'].nunique()
            counts = counts.fillna(0).astype(int)
            counts['donors prop.'] = counts['donors'] / ndonors
            counts['donors prop.'] = counts['donors prop.'].apply(lambda x: f"{x*100:.1f} %")
            counts = counts.sort_values(by=['donors', 'total'], ascending=[False, False])
            print(counts.head(20))

        print()
        print('------------------------------------------')
        print('--- Mutation Enrichment (top 20 genes) ---')
        print('------------------------------------------')
        print()
        print(rframe.set_index('gene').head(20))

    sys.stdout = original_stdout


def load_cmdline_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Simple analysis and summaries for single genes.')
    parser.add_argument('--posmuts', type=str, required=True, help='Path to positive class mutations (tsv).')
    parser.add_argument('--negmuts', type=str, required=True, help='Path to negative class mutations (tsv).')
    parser.add_argument('--outfile', type=str, required=True, help='Output file basename (no extension pls).')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    main()