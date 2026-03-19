
import os
import sys
import argparse
import numpy as np
import pandas as pd 
from statsmodels.stats.multitest import multipletests
from statsmodels.formula.api import logit
from glob import glob 

import statsmodels.formula.api as smf

from fileio import load_mutations
from fileio import load_gmt
from selection import select_genesets

from reporting import (
    summarise_basic_info,
    summarise_vclasses,
    summarise_donors,
)

from filtering import (
    filter_donors_without_prostate,
    filter_donors_without_dpclust,
    filter_donors_without_trees,
    filter_hypermutators
)

from formatting import generate_geneset_matrix

import warnings
warnings.filterwarnings('ignore')

pd.options.display.float_format = "{:,.3f}".format
pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 100)
pd.set_option('display.width', 150)

INFILE_MUTATIONS = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/manual/mutations.assigned.160326.tsv'
INFILE_SHEET = '/home/grace/work/PPCG_DifferentialGenesetMutation/samplesheet.angel.alldonors.tsv'
INFILE_GENESETS = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/genesets/c2.cp.kegg_medicus.v2026.1.Hs.symbols.gmt'
# INFILE_GENESETS = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/genesets/c2.cp.reactome.v2026.1.Hs.symbols.gmt'
# INFILE_GENESETS = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/genesets/h.all.v2026.1.Hs.symbols.gmt'
# INFILE_GENESETS = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/genesets/c5.go.bp.v2026.1.Hs.symbols.gmt'
INDIR_DPCLUST = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/angel/dpclust_ccf'
INDIR_CONIPHER = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/angel/conipher_trees'
INDIR_TIMING = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/angel/time_trees_acceleration_5'

def main():
    args = load_cmdline_args()
    muts = load_mutations(args.muts, args.sheet)
    muts = filter_mutations(muts, args)
    report_mutations(muts, args)

    # geneset selection 
    gset_LUT = load_gmt(args.gmt)
    gset_LUT, _ = select_genesets(gset_LUT, muts)

    # mutation matrix 
    mat = generate_geneset_matrix(gset_LUT, muts, args)
    mat = mat.drop('metastatic', axis=1)

    # add burden
    burden_LUT = muts.groupby('donor')['gene'].nunique().to_dict()
    mat['burden'] = burden_LUT

    # add latency
    latency_LUT = load_latency(sorted(list(muts['donor'].unique())), args)
    mat['latency'] = latency_LUT

    # reformat column order
    prepend = ['burden', 'latency']
    mat = mat[prepend+[x for x in mat.columns if x not in prepend]].copy()

    # remove donors with only a few trunk mutations
    mat = mat[mat['burden']>=10].copy()
    print(mat.shape)
    
    res = analyze_gene_associations(mat)
    print()
    print(res.head(10))
    
    query = res['gene'].values[0]
    print()
    print(mat[prepend+[query]].sort_values('latency'))



def analyze_gene_associations(df: pd.DataFrame, alpha: float = 0.05) -> pd.DataFrame:
    """
    Test the association between each gene mutation and a continuous response variable,
    adjusting for mutational burden.

    The model fit for each gene is:
        response ~ intercept + mutation_gene + burden

    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe where:
            - Column 0 : 'burden'   (int)   total number of mutations per patient
            - Column 1 : 'response' (float) continuous outcome variable
            - Columns 2+: gene columns (0 = no mutation, 1 = mutation)
    alpha : float
        Significance threshold for the 'significant' flag (default 0.05).

    Returns
    -------
    pd.DataFrame
        One row per gene with columns:
            gene         : gene name
            n_mutated    : number of patients with a mutation in this gene
            coef         : regression coefficient for mutation (burden-adjusted effect)
            se           : standard error of the coefficient
            ci_lower     : lower bound of 95% confidence interval
            ci_upper     : upper bound of 95% confidence interval
            t_stat       : t-statistic
            p_value      : two-sided p-value
            p_value_fdr  : Benjamini-Hochberg FDR-adjusted p-value
            significant  : bool — True if p_value_fdr < alpha
            direction    : 'higher', 'lower', or 'n.s.'
        Sorted by p_value ascending.
    """
    burden_col = df.columns[0]
    response_col = df.columns[1]
    gene_cols = df.columns[2:]

    rows = []

    for gene in gene_cols:
        n_mutated = df[gene].sum()

        # Skip genes with no variation (all 0s or all 1s)
        if n_mutated == 0 or n_mutated == len(df):
            rows.append({
                "gene": gene,
                "n_mutated": int(n_mutated),
                "coef": np.nan, "se": np.nan,
                "ci_lower": np.nan, "ci_upper": np.nan,
                "t_stat": np.nan, "p_value": np.nan,
                "p_value_fdr": np.nan, "significant": False,
                "direction": "n.s. (no variance)",
            })
            continue

        model = smf.ols(
            formula=f"response ~ mutation + burden",
            data=df.rename(columns={gene: "mutation", burden_col: "burden", response_col: "response"}),
        ).fit()

        coef = model.params["mutation"]
        se = model.bse["mutation"]
        t_stat = model.tvalues["mutation"]
        p_value = model.pvalues["mutation"]
        ci_lower, ci_upper = model.conf_int().loc["mutation"]

        rows.append({
            "gene": gene,
            "n_mutated": int(n_mutated),
            "coef": round(coef, 4),
            "se": round(se, 4),
            "ci_lower": round(ci_lower, 4),
            "ci_upper": round(ci_upper, 4),
            "t_stat": round(t_stat, 4),
            "p_value": round(p_value, 6),
            "p_value_fdr": np.nan,  # filled below
            "significant": False,
            "direction": "",
        })

    results = pd.DataFrame(rows)

    # --- Benjamini-Hochberg FDR correction on testable genes only ---
    testable = results["p_value"].notna()
    if testable.sum() > 0:
        p_vals = results.loc[testable, "p_value"].values
        fdr = _bh_correction(p_vals)
        results.loc[testable, "p_value_fdr"] = np.round(fdr, 6)
        results.loc[testable, "significant"] = fdr < alpha
        results["direction"] = results.apply(
            lambda r: ("higher" if r["coef"] > 0 else "lower") if r["significant"] else "n.s.",
            axis=1,
        )

    return results.sort_values("p_value").reset_index(drop=True)


def _bh_correction(p_values: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction. Returns adjusted p-values."""
    n = len(p_values)
    order = np.argsort(p_values)
    ranks = np.argsort(order) + 1          # 1-based ranks
    adjusted = p_values * n / ranks
    # Enforce monotonicity (step-down)
    adjusted_ordered = adjusted[order]
    for i in range(n - 2, -1, -1):
        adjusted_ordered[i] = min(adjusted_ordered[i], adjusted_ordered[i + 1])
    fdr = np.empty(n)
    fdr[order] = adjusted_ordered
    return np.minimum(fdr, 1.0)



def select_genes(muts: pd.DataFrame) -> set[str]:
    TOP_N = 50
    counts = muts[muts['cohort']=='COMBI'].groupby('gene')['donor'].nunique()
    counts = counts.sort_values(ascending=False)
    print(counts.head())
    selected = set(counts.head(TOP_N).index.to_list())
    return selected

def load_latency(donors: list[str], args: argparse.Namespace) -> dict[str, float]:
    out = {}
    for donor in donors:
        filepath = f"{args.timing_dir}/{donor}_timing_snv_model_acceleration.tsv"
        df = pd.read_csv(filepath, sep='\t', header=0)
        earliest_clone = df['model_median_time'].min()
        latest_clone = df['model_median_time'].max()
        out[donor] = latest_clone-earliest_clone
    return out

def filter_mutations(df: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    print('filtering donors/mutations')
    df = df[df['cohort']=='COMBI'].copy()
    # df = df[df['vclass']!='CNA'].copy()
    print(df['donor'].nunique(), df['ID'].nunique())
    
    df = filter_donors_without_prostate(df)
    # df = filter_donors_without_dpclust(df, args.dpclust_dir)
    df = filter_donors_without_trees(df, args.conipher_dir)

    # has timing data
    TIMING_DONORS = set([x.split('/')[-1][:8] for x in glob(f"{INDIR_TIMING}/*.tsv")])
    df = df[df['donor'].isin(TIMING_DONORS)]
    print(df['donor'].nunique(), df['ID'].nunique())

    # only keep trunk mutations
    df = df[df['asmt']=='Primary:Trunk'].copy()
    print(df['donor'].nunique(), df['ID'].nunique())

    df = filter_hypermutators(df, max_mutated_genes=args.hypermutator_ngenes)
    print(df['donor'].nunique(), df['ID'].nunique())
    return df 


def report_mutations(muts: pd.DataFrame, args: argparse.Namespace) -> None:
    # outfile_report = f"{args.outdir}/{args.run_id}/mutation_summary.txt"
    
    # print(f'writing report to {outfile_report}')
    # original_stdout = sys.stdout
    # with open(outfile_report, 'w') as f:
    #     sys.stdout = f
    summarise_basic_info(muts)
    summarise_vclasses(muts)
    summarise_donors(muts)
    # sys.stdout = original_stdout


def load_cmdline_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Runs downstream mutation analysis.')
    
    # run ID
    parser.add_argument('--run-id', type=str, default='default', help='ID for this analysis.')
    parser.add_argument('--outdir', type=str, default='./results', help='Path to output directory.')

    # files
    parser.add_argument('--gmt', type=str, default=INFILE_GENESETS, help='Path to msigdb genesets .gmt file.')
    parser.add_argument('--muts', type=str, default=INFILE_MUTATIONS, help='Path to assigned mutations file.')
    parser.add_argument('--sheet', type=str, default=INFILE_SHEET, help='Path to samplesheet (metadata).')
    parser.add_argument('--timing-dir', type=str, default=INDIR_TIMING, help='Path to timing folder')
    parser.add_argument('--dpclust-dir', type=str, default=INDIR_DPCLUST, help='Path to DPClust folder')
    parser.add_argument('--conipher-dir', type=str, default=INDIR_CONIPHER, help='Path to conipher folder')
    
    # donor selection
    parser.add_argument('--hypermutator-ngenes', type=int, default=600, help="Number of mutated genes defining a hypermutator")

    args = parser.parse_args()
    return args


if __name__ == '__main__':
    main()
