
import os
import sys
import argparse
import numpy as np
import pandas as pd 
from scipy.stats import fisher_exact, chi2
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
pd.set_option('display.width', None)

INFILE_MUTATIONS = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/manual/mutations.assigned.160326.tsv'
INFILE_SHEET = '/home/grace/work/PPCG_DifferentialGenesetMutation/samplesheet.angel.alldonors.tsv'
# INFILE_GENESETS = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/genesets/c2.cp.kegg_medicus.v2026.1.Hs.symbols.gmt'
# INFILE_GENESETS = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/genesets/c2.cp.reactome.v2026.1.Hs.symbols.gmt'
INFILE_GENESETS = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/genesets/h.all.v2026.1.Hs.symbols.gmt'
# INFILE_GENESETS = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/genesets/c5.go.bp.v2026.1.Hs.symbols.gmt'
INDIR_DPCLUST = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/angel/dpclust_ccf'
INDIR_CONIPHER = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/angel/conipher_trees'
INDIR_TIMING = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/angel/time_trees_acceleration_5'

def main():
    args = load_cmdline_args()
    rundir = f"{args.outdir}/{args.run_id}"
    if not os.path.exists(rundir):
        os.makedirs(rundir, exist_ok=True)

    ### mutations
    muts = load_mutations(args.muts, args.sheet)
    muts = filter_mutations(muts, args)
    report_mutations(muts, args)

    ### covariate table
    all_donors = sorted(list(muts['donor'].unique()))
    lframe = pd.DataFrame(index=all_donors)

    # add burden
    burden_LUT = muts.groupby('donor')['gene'].nunique().to_dict()
    lframe['burden'] = burden_LUT

    # add latency
    latency_LUT = load_latency(sorted(list(muts['donor'].unique())), args)
    lframe['latency'] = latency_LUT
    lframe['log_latency'] = np.log(lframe['latency'])

    # add groups
    lframe = lframe.sort_values('latency')
    mid = lframe.shape[0]//2
    group_LUT = {}
    for i, donor in enumerate(lframe.index.to_list()):
        group_LUT[donor] = 'FAST' if i <= mid else 'SLOW'
    lframe['group'] = group_LUT 

    ### mutation matrix 
    if args.genes:
        genes = select_genes(muts)
        gset_LUT = {gene: [gene] for gene in genes}
    else:
        gset_LUT = load_gmt(args.gmt)
        gset_LUT, sframe = select_genesets(gset_LUT, muts)

    mat = generate_geneset_matrix(gset_LUT, muts, all_donors, args)
    # mat = pd.concat([lframe[['burden', 'latency', 'log_latency', 'group']], mat], axis=1)

    ### continuous assessment
    df = pd.concat([lframe[['burden', 'log_latency']], mat], axis=1).copy()
    res_cont = analyze_continuous(df)
    res_cont = res_cont.rename(columns={'gene': 'geneset'})
    
    ### categorical assessment
    # format matrix
    df = pd.concat([lframe[['burden', 'group']], mat], axis=1).copy()
    df['group'] = df['group'].map({'FAST': 1, 'SLOW': 0}).astype(int)
    
    # logistic test
    res_logit = analyze_logistic(df, label_col='group')
    
    # fisher test
    res_fish = analyze_fisher(df, label_col='group')
    res_fish = res_fish.rename(columns={'gene': 'geneset'})

    if args.genes:
        report_results(res_cont, res_logit, res_fish, lframe, mat, args, sframe=None)
    else:
        report_results(res_cont, res_logit, res_fish, lframe, mat, args, sframe=sframe)

    # import matplotlib.pyplot as plt 
    # df['log_latency'] = np.log(df['latency'])
    # plt.scatter(x=list(range(df.shape[0])), y=df['latency'].to_list())
    # plt.savefig('/home/grace/work/PPCG_DifferentialGenesetMutation/results/latency_continuous/test/latency.png')
    # plt.close()
    # plt.scatter(x=list(range(df.shape[0])), y=df['log_latency'].to_list())
    # plt.savefig('/home/grace/work/PPCG_DifferentialGenesetMutation/results/latency_continuous/test/log_latency.png')
    # plt.close()

def select_genes(muts: pd.DataFrame) -> set[str]:
    TOP_N = 20
    counts = muts.groupby('gene')['donor'].nunique()
    counts = counts.sort_values(ascending=False)
    selected = set(counts.head(TOP_N).index.to_list())
    return selected

def report_results(
    res_cont: pd.DataFrame, 
    res_logit: pd.DataFrame, 
    res_fish: pd.DataFrame, 
    lframe: pd.DataFrame, 
    mat: pd.DataFrame, 
    args: argparse.Namespace,
    sframe: pd.DataFrame|None, 
    ) -> None:
    
    outfile_continuous = f"{args.outdir}/{args.run_id}/results_continuous.tsv"
    outfile_fisher = f"{args.outdir}/{args.run_id}/results_fisher.tsv"
    outfile_logistic = f"{args.outdir}/{args.run_id}/results_logistic.tsv"
    outfile_summary = f"{args.outdir}/{args.run_id}/results_summary.txt"

    items_cont = set(res_cont['geneset'].unique())
    items_fish = set(res_fish['geneset'].unique())
    items_logit = set(res_logit['geneset'].unique())
    assert items_cont == items_fish
    assert items_cont == items_logit
    assert items_fish == items_logit

    res_cont = res_cont.set_index('geneset')
    df = pd.concat([lframe, mat], axis=1).copy()
    for gset, row in res_cont.iterrows():
        res_cont.loc[gset, 'med_latency_mut'] = df[df[gset]==1]['latency'].median()
        res_cont.loc[gset, 'med_latency_non_mut'] = df[df[gset]==0]['latency'].median()
    res_cont = res_cont.reset_index()

    res_cont = res_cont.sort_values('p_value')
    res_logit = res_logit.sort_values('p_value')
    res_fish = res_fish.sort_values('p_value')

    res_fish = res_fish.rename(columns={'n_mutated_meta': 'n_mutated_FAST', 'n_mutated_nonmeta': 'n_mutated_SLOW'})
    res_logit = res_logit.rename(columns={'n_mutated_meta': 'n_mutated_FAST', 'n_mutated_nonmeta': 'n_mutated_SLOW'})
    
    res_cont.to_csv(outfile_continuous, sep='\t', index=False, float_format='%.5f')
    res_logit.to_csv(outfile_logistic, sep='\t', index=False, float_format='%.5f')
    res_fish.to_csv(outfile_fisher, sep='\t', index=False, float_format='%.5f')

    original_stdout = sys.stdout
    with open(outfile_summary, 'w') as f:
        sys.stdout = f

        print()
        print(f"DONORS: {lframe.shape[0]}")
        print(f"GENESETS: {len(items_cont)}")

        print()
        if sframe is not None:
            tophits = set() 
            tophits = tophits | set(res_cont.head(10)['geneset'].unique()) 
            tophits = tophits | set(res_fish.head(10)['geneset'].unique()) 
            tophits = tophits | set(res_logit.head(10)['geneset'].unique())
            print()
            print('[Geneset Selection]')
            print(sframe[sframe['geneset'].isin(tophits)].set_index('geneset'))

        print()
        print('--------------')
        print('--- FISHER ---')
        print('--------------')
        print(res_fish.set_index('geneset').head(10))
        
        print()
        print('----------------')
        print('--- LOGISTIC ---')
        print('----------------')
        print(res_logit.set_index('geneset').head(10))
        
        print()
        print('------------------')
        print('--- CONTINUOUS ---')
        print('------------------')

        fields = ['n_mutated', 'med_latency_mut', 'med_latency_non_mut', 'p_value', 'p_value_fdr', 'direction']
        print(res_cont.set_index('geneset')[fields].head(10))

        # df = mat[['burden', 'latency']].sort_values('latency').copy()
        df = lframe.copy()
        genesets = res_cont.head(10)['geneset'].to_list()
        for i, gset in enumerate(genesets):
            df[i+1] = mat[gset]
        print()
        print(df) 
        print()
        for i, gset in enumerate(genesets):
            print(f"{i+1}: {gset}")
        print()

    sys.stdout = original_stdout




def analyze_fisher(df: pd.DataFrame,
                                burden_col: str = "burden",
                                label_col: str = "label",
                                fdr_threshold: float = 0.05) -> pd.DataFrame:
    """
    Test each gene for association with a binary outcome using Fisher's exact test,
    with Benjamini-Hochberg FDR correction for multiple comparisons.

    Parameters
    ----------
    df : pd.DataFrame
        One row per patient. label_col is binary (0/1); all other columns are
        binary gene mutation status (0/1).
    label_col : str
        Name of the outcome column.
    fdr_threshold : float
        FDR threshold for the gene-level correction.

    Returns
    -------
    pd.DataFrame
        Columns: gene, n_mutated_meta, n_mutated_nonmeta, odds_ratio,
                 p_value, q_value, significant
    """
    gene_cols = [c for c in df.columns if c != label_col and c != burden_col]
    meta      = df[label_col] == 1
    non_meta  = df[label_col] == 0

    results = []
    for gene in gene_cols:
        mut = df[gene] == 1
        a = (meta     &  mut).sum()
        b = (meta     & ~mut).sum()
        c = (non_meta &  mut).sum()
        d = (non_meta & ~mut).sum()
        odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative="two-sided")
        results.append({
            "gene":              gene,
            "n_mutated_meta":    int(a),
            "n_mutated_nonmeta": int(c),
            "odds_ratio":        odds_ratio,
            "p_value":           p_value,
        })

    results_df = pd.DataFrame(results)
    reject, q_values, _, _ = multipletests(
        results_df["p_value"], method="fdr_bh", alpha=fdr_threshold
    )
    results_df["q_value"]     = q_values
    results_df["significant"] = reject
    return results_df.sort_values("odds_ratio", ascending=False).reset_index(drop=True)




def analyze_logistic(
    df: pd.DataFrame, 
    burden_col: str = "burden", 
    label_col: str = "label", 
    fdr_threshold: float = 0.05) -> pd.DataFrame:
    """
    Test each geneset for association with a binary outcome using logistic regression,
    adjusting for per-patient mutational burden.

    Model per geneset:
        label ~ geneset_mutated + burden

    Parameters
    ----------
    df : pd.DataFrame
        One row per patient. label_col is binary (0/1); all other columns are
        binary geneset mutation status (0/1).
    label_col : str
        Name of the outcome column.
    fdr_threshold : float
        FDR threshold for the geneset-level correction.

    Returns
    -------
    pd.DataFrame
        Columns: geneset, n_mutated_meta, n_mutated_nonmeta, log_odds_ratio,
                 odds_ratio, p_value, q_value, significant
    """
    gset_cols = [c for c in df.columns if c != label_col and c != burden_col]
    total_muts = df[burden_col]
    meta       = df[label_col] == 1
    non_meta   = df[label_col] == 0

    results = []
    for gset in gset_cols:
        burden = total_muts
        burden_std = burden.std()
        burden_scaled = (burden - burden.mean()) / burden_std if burden_std > 0 else burden - burden.mean()
        mut = df[gset]
        a = (meta     & (mut == 1)).sum()
        c = (non_meta & (mut == 1)).sum()

        if mut.nunique() < 2:
            results.append({
                "geneset": gset, "n_mutated_meta": int(a), "n_mutated_nonmeta": int(c),
                "log_odds_ratio": np.nan, "odds_ratio": np.nan, "p_value": np.nan,
            })
            continue

        fit_df = pd.DataFrame({
            label_col:  df[label_col].values,
            "gset_mut": mut.values,
            "burden":   burden_scaled.values,
        })
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                model   = logit(f"{label_col} ~ gset_mut + burden", data=fit_df).fit(
                    disp=False, maxiter=200)
            log_or  = model.params["gset_mut"]
            p_value = model.pvalues["gset_mut"]
        except Exception:
            log_or, p_value = np.nan, np.nan

        results.append({
            "geneset":           gset,
            "n_mutated_meta":    int(a),
            "n_mutated_nonmeta": int(c),
            "log_odds_ratio":    log_or,
            "odds_ratio":        np.exp(log_or) if not np.isnan(log_or) else np.nan,
            "p_value":           p_value,
        })

    results_df = pd.DataFrame(results)
    valid   = results_df["p_value"].notna()
    reject  = np.zeros(len(results_df), dtype=bool)
    q_vals  = np.full(len(results_df), np.nan)
    if valid.sum() > 0:
        rej_v, q_v, _, _ = multipletests(
            results_df.loc[valid, "p_value"], method="fdr_bh", alpha=fdr_threshold)
        reject[valid] = rej_v
        q_vals[valid] = q_v
    results_df["q_value"]     = q_vals
    results_df["significant"] = reject
    return results_df.sort_values("odds_ratio", ascending=False).reset_index(drop=True)

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

    # hypermutators
    df = filter_hypermutators(df, max_mutated_genes=args.hypermutator_ngenes)
    print(df['donor'].nunique(), df['ID'].nunique())

    counts = df.groupby('donor')['gene'].nunique()
    valid = set(counts[counts>=10].index.to_list())
    df = df[df['donor'].isin(valid)]
    print(df['donor'].nunique(), df['ID'].nunique())

    return df 

def report_mutations(muts: pd.DataFrame, args: argparse.Namespace) -> None:
    outfile_report = f"{args.outdir}/{args.run_id}/mutations_summary.txt"
    
    print(f'writing report to {outfile_report}')
    original_stdout = sys.stdout
    with open(outfile_report, 'w') as f:
        sys.stdout = f
        summarise_basic_info(muts)
        summarise_vclasses(muts)
        summarise_donors(muts)
        summarise_assignments(muts)
    sys.stdout = original_stdout

def summarise_assignments(muts: pd.DataFrame) -> None:
    print('\n--- Mutation Assignments ---\n')
    print()
    print(muts.drop_duplicates('ID')['asmt'].value_counts())
    print()
    print(muts.drop_duplicates('ID').groupby('donor')['asmt'].value_counts())

def load_latency(donors: list[str], args: argparse.Namespace) -> dict[str, float]:
    out = {}
    for donor in donors:
        filepath = f"{args.timing_dir}/{donor}_timing_snv_model_acceleration.tsv"
        df = pd.read_csv(filepath, sep='\t', header=0)
        earliest_clone = df['model_median_time'].min()
        latest_clone = df['model_median_time'].max()
        out[donor] = latest_clone-earliest_clone
    return out

def analyze_continuous(df: pd.DataFrame, alpha: float = 0.05) -> pd.DataFrame:
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

    # genes not genesets
    parser.add_argument('--genes', action="store_true", help="Whether to run genes instead of genesets")

    args = parser.parse_args()
    return args


if __name__ == '__main__':
    main()
