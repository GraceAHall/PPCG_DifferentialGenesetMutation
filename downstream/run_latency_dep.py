
import os
import sys
import argparse
import numpy as np
import pandas as pd 
from collections import defaultdict
from statsmodels.stats.multitest import multipletests
from statsmodels.formula.api import logit
from glob import glob 
from fileio import load_mutations
from fileio import load_gmt

from filtering import (
    filter_donors_without_prostate,
    filter_donors_without_dpclust,
    filter_donors_without_trees,
    filter_hypermutators
)
from reporting import (
    summarise_basic_info,
    summarise_vclasses,
    summarise_donors,
)
# from selection import select_genesets
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
# INFILE_GENESETS = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/genesets/c5.go.bp.v2026.1.Hs.symbols.gmt'
# INFILE_GENESETS = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/genesets/h.all.v2026.1.Hs.symbols.gmt'
# INFILE_GENESETS = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/genesets/c2.cp.reactome.v2026.1.Hs.symbols.gmt'
# INFILE_GENESETS = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/genesets/c5.go.bp.v2026.1.Hs.symbols.gmt'
INFILE_GENESETS = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/genesets/c2.cp.kegg_medicus.v2026.1.Hs.symbols.gmt'
INDIR_DPCLUST = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/angel/dpclust_ccf'
INDIR_CONIPHER = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/angel/conipher_trees'
INDIR_TIMING = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/angel/time_trees_acceleration_5'

def main():
    args = load_cmdline_args()
    muts = load_mutations(args.muts, args.sheet)
    # muts = muts[muts['vclass']!='CNA']
    muts = filter_mutations(muts, args)
    report_mutations(muts, args)

    gset_LUT = load_gmt(args.gmt)
    # gset_LUT, sel_df = select_genesets(gset_LUT, muts)
    gset_LUT = select_genesets(gset_LUT, muts)
    all_genesets = sorted(list(gset_LUT.keys()))
    burden_LUT = muts.groupby('donor')['gene'].nunique().to_dict()
    latency_LUT = load_latency(sorted(list(muts['donor'].unique())), args)
    
    print(f"\nTotal genesets to test: {len(gset_LUT)}")

    mat = generate_geneset_matrix(gset_LUT, muts, args)
    mat = mat.drop('metastatic', axis=1)
    mat['burden'] = burden_LUT
    burden = mat['burden']
    burden_std = burden.std()
    burden_scaled = (burden - burden.mean()) / burden_std if burden_std > 0 else burden - burden.mean()
    mat['burden_scaled'] = burden_scaled

    mat['latency'] = latency_LUT
    latency = mat['latency']
    latency_std = latency.std()
    latency_scaled = (latency - latency.mean()) / latency_std if latency_std > 0 else latency - latency.mean()
    mat['latency_scaled'] = latency_scaled

    mat = mat[['burden', 'burden_scaled', 'latency', 'latency_scaled']+all_genesets].copy()
    mat = mat[mat['burden']>=10].copy()
    res = run_logit_gset_latency_association(mat, all_genesets)
    res = res.sort_values('p_value')
    print(res['significant'].value_counts())
    print()
    print(res.head(10))

    gset = res['geneset'].values[0]
    genes = set(gset_LUT[gset])
    df = mat[['burden', 'latency']].sort_values('latency').copy()
    df['Trunk'] = mat[gset]
    
    muts = load_mutations(args.muts, args.sheet)
    muts = muts[muts['donor'].isin(set(df.index.to_list()))].copy()
    print(muts['asmt'].value_counts())
    muts = muts[~muts['asmt'].isin(['Primary:Trunk'])].copy()
    # muts = muts[~muts['asmt'].isin(['Primary:Leaf', 'Primary:Trunk', 'Primary:Branch'])].copy()
    muts = muts[muts['gene'].isin(genes)].copy()
    mut_donors = set(muts['donor'].unique())
    for donor in df.index.to_list():
        df.loc[donor, 'Metastatic'] = 1 if donor in mut_donors else 0
    df['Metastatic'] = df['Metastatic'].astype(int)
    print()
    print(gset)
    print(df)
    print()

def run_logit_gset_latency_association(
    df: pd.DataFrame,
    all_genesets: list[str],
    fdr_threshold: float = 0.05) -> pd.DataFrame:
    """
    Test each geneset for association with a binary outcome using logistic regression,
    adjusting for per-patient mutational burden.

    Model per geneset:
        metastatic ~ geneset_mutated + burden

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
    gset_cols  = all_genesets

    results = []
    for gset in gset_cols:
        burden_scaled = df['burden_scaled']
        latency_scaled = df['latency_scaled']
        mut = df[gset]
        ndonors = mut.sum()

        fit_df = pd.DataFrame({
            "gset_mut": mut.values,
            "latency":  latency_scaled.values,
            "burden":   burden_scaled.values,
        })
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                model   = logit(f"gset_mut ~ latency + burden", data=fit_df).fit(
                    disp=False, maxiter=200)
            log_or  = model.params["latency"]
            p_value = model.pvalues["latency"]
        except Exception:
            log_or, p_value = np.nan, np.nan

        results.append({
            "geneset":           gset,
            "mut_donors":        ndonors,
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

def load_latency(donors: list[str], args: argparse.Namespace) -> dict[str, float]:
    out = {}
    for donor in donors:
        filepath = f"{args.timing_dir}/{donor}_timing_snv_model_acceleration.tsv"
        df = pd.read_csv(filepath, sep='\t', header=0)
        earliest_clone = df['model_median_time'].min()
        latest_clone = df['model_median_time'].max()
        out[donor] = latest_clone-earliest_clone
    return out

def select_genesets(gset_LUT: dict[str, list[str]], muts: pd.DataFrame) -> dict[str, list[str]]:
    MIN_OVERLAP = 5
    MIN_DONORS = 10

    # df = muts[muts['cohort']=='COMBI'].copy()
    df = muts.copy()
    ignore = set()
    for gset, genes in gset_LUT.items():
        dfslice = df[df['gene'].isin(set(genes))]
        ngenes = dfslice['gene'].nunique()
        ndonors = dfslice['donor'].nunique()
        if ngenes < MIN_OVERLAP or ndonors < MIN_DONORS:
            ignore.add(gset)

    gset_LUT = {k: v for k, v in gset_LUT.items() if k not in ignore}
    return gset_LUT


def filter_mutations(df: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    print('filtering donors/mutations')
    df = df[df['cohort']=='COMBI'].copy()
    
    # df = filter_donors_without_prostate(df)
    # df = filter_donors_without_dpclust(df, args.dpclust_dir)
    # df = filter_donors_without_trees(df, args.conipher_dir)

    # has timing data
    TIMING_DONORS = set([x.split('/')[-1][:8] for x in glob(f"{INDIR_TIMING}/*.tsv")])
    df = df[df['donor'].isin(TIMING_DONORS)]

    # only keep trunk mutations
    print(df['donor'].nunique(), df['ID'].nunique())
    df = df[df['asmt']=='Primary:Trunk'].copy()
    # df = df[df['asmt'].isin(['Primary:Trunk', 'Primary:Branch'])].copy()
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
    parser.add_argument('--muts', type=str, default=INFILE_MUTATIONS, help='Path to assigned mutations file.')
    parser.add_argument('--sheet', type=str, default=INFILE_SHEET, help='Path to samplesheet (metadata).')
    parser.add_argument('--gmt', type=str, default=INFILE_GENESETS, help='Path to msigdb genesets .gmt file.')
    parser.add_argument('--timing-dir', type=str, default=INDIR_TIMING, help='Path to timing folder')
    parser.add_argument('--dpclust-dir', type=str, default=INDIR_DPCLUST, help='Path to DPClust folder')
    parser.add_argument('--conipher-dir', type=str, default=INDIR_CONIPHER, help='Path to conipher folder')
    
    # donor selection
    parser.add_argument('--hypermutator-ngenes', type=int, default=500, help="Number of mutated genes defining a hypermutator")

    args = parser.parse_args()
    return args


if __name__ == '__main__':
    main()
