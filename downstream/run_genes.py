
import os
import sys
import argparse
import numpy as np
import pandas as pd 
from collections import defaultdict
from statsmodels.stats.multitest import multipletests
from scipy.stats import fisher_exact

from fileio import load_mutations
from fileio import load_gmt

from filtering import (
    filter_ppcg_control,
    filter_donors_without_prostate,
    filter_donors_without_dpclust,
    filter_donors_without_trees,
    filter_primary_leaves,
    filter_secondary_leaves,
    filter_hypermutators
)

from reporting import (
    summarise_basic_info,
    summarise_vclasses,
    summarise_donors,
    plot_distribution
)

from formatting import generate_gene_matrix
from stats import run_fisher_gene_association
from stats import run_logit_gene_association

import warnings
warnings.filterwarnings('ignore')

pd.options.display.float_format = "{:,.3f}".format
pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 200)
pd.set_option('display.width', 180)

INFILE_MUTATIONS = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/manual/mutations.assigned.160326.tsv'
INFILE_SHEET = '/home/grace/work/PPCG_DifferentialGenesetMutation/samplesheet.angel.alldonors.tsv'
INDIR_DPCLUST = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/angel/dpclust_ccf'
INDIR_CONIPHER = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/angel/conipher_trees'


def main():
    args = load_cmdline_args()
    outdir = f"{args.outdir}/{args.run_id}"
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    report_settings(args)

    # mutations
    muts = load_mutations(args.muts, args.sheet)

    # filtering
    muts = filter_mutations_alt(muts, args)

    # reporting
    report_mutations(muts, args)

    # gene selection 
    genes = select_genes(muts)

    # mutation matrix 
    mat = generate_gene_matrix(muts, genes)
    print(mat.shape)
    print(mat.iloc[:10, :6])
    print()
    print(mat.iloc[:, 2:].sum().sort_values(ascending=False).head())
    
    # statistical tests
    do_gene_enrichment(mat, args)

def do_gene_enrichment(df: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    print('running gene enrichment tests')
    outfile_fish = f"{args.outdir}/{args.run_id}/gene_fisher.tsv"
    outfile_logit = f"{args.outdir}/{args.run_id}/gene_logit.tsv"
    outfile_comb = f"{args.outdir}/{args.run_id}/gene_combined.tsv"
    
    res_fish = run_fisher_gene_association(df)
    res_logit = run_logit_gene_association(df)

    # res_fish['odds_ratio'] = res_fish['odds_ratio'].clip(upper=100)
    # res_logit['odds_ratio'] = res_logit['odds_ratio'].clip(upper=100)
    # res_logit['log_odds_ratio'] = res_logit['log_odds_ratio'].clip(upper=100)
    res_fish = res_fish.sort_values('p_value')
    res_logit = res_logit.sort_values('p_value')

    res_fish.to_csv(outfile_fish, sep='\t', index=False, float_format='%.5f')
    res_logit.to_csv(outfile_logit, sep='\t', index=False, float_format='%.5f')

    res_fish = res_fish.set_index('gene')
    res_logit = res_logit.set_index('gene')
    res_comb = res_logit.copy()
    mask = (res_comb['n_mutated_meta']==0) | (res_comb['n_mutated_nonmeta']==0)
    res_comb.loc[mask, 'odds_ratio'] = res_fish.loc[mask, 'odds_ratio']
    res_comb.loc[mask, 'p_value'] = res_fish.loc[mask, 'p_value']
    res_comb.loc[mask, 'q_value'] = res_fish.loc[mask, 'q_value']
    res_comb.loc[mask, 'significant'] = res_fish.loc[mask, 'significant']
    res_comb = res_comb.reset_index()
    res_comb = res_comb.sort_values('p_value')
    
    res_comb.to_csv(outfile_comb, sep='\t', index=False, float_format='%.5f')
    return res_comb

def select_genes(muts: pd.DataFrame) -> set[str]:
    TOP_N = 50
    counts = muts[muts['cohort']=='COMBI'].groupby('gene')['donor'].nunique()
    counts = counts.sort_values(ascending=False)
    print(counts.head())
    selected = set(counts.head(TOP_N).index.to_list())
    return selected

def report_settings(args: argparse.Namespace) -> None:
    outfile_settings = f"{args.outdir}/{args.run_id}/settings.txt"
    
    original_stdout = sys.stdout
    with open(outfile_settings, 'w') as f:
        sys.stdout = f
        print('Arguments for this analysis ---\n')
        for arg_name, arg_value in vars(args).items():
            print(f"- {arg_name}: {arg_value}")
    sys.stdout = original_stdout

def report_mutations(muts: pd.DataFrame, args: argparse.Namespace) -> None:
    df = muts.copy()
    outfile_report = f"{args.outdir}/{args.run_id}/mutation_summary.txt"
    outfile_plot = f"{args.outdir}/{args.run_id}/mutation_hist.png"
    
    print(f'writing report to {outfile_report}')
    original_stdout = sys.stdout
    with open(outfile_report, 'w') as f:
        sys.stdout = f
        summarise_basic_info(df)
        summarise_vclasses(df)
        summarise_donors(df)
    sys.stdout = original_stdout

    print(f'writing plot to {outfile_plot}')
    plot_distribution(df, filepath=outfile_plot)

def filter_mutations_alt(df: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    print('filtering donors/mutations')
    df = df[df['vclass']!='CNA'].copy()
    df_ppcg = df[df['cohort']=='PPCG'].copy()
    df_combi = df[df['cohort']=='COMBI'].copy()
    
    df_combi = filter_donors_without_prostate(df_combi)
    # df_combi = filter_donors_without_dpclust(df_combi, args.dpclust_dir)
    df_combi = filter_donors_without_trees(df_combi, args.conipher_dir)
    df_combi = df_combi[df_combi['asmt'].isin(['Primary:Trunk', 'Primary:Branch', 'Primary:Seed'])].copy()
    df_combi = filter_hypermutators(df_combi, max_mutated_genes=args.hypermutator_ngenes)
    df_ppcg = filter_hypermutators(df_ppcg, max_mutated_genes=args.hypermutator_ngenes)
    df_ppcg = pick_matched_cohort(df_ppcg, df_combi)

    # if not args.ppcg_all:
        # df_ppcg = filter_ppcg_control(df_ppcg)
    # if args.filter_primary_leaves:
    #     df_combi = filter_primary_leaves(df_combi, verbose=True)
    # if args.filter_secondary_leaves:
    #     df_combi = filter_secondary_leaves(df_combi, verbose=True)

    df = pd.concat([df_ppcg, df_combi], ignore_index=True)
    return df 

def pick_matched_cohort(df_ppcg: pd.DataFrame, df_combi: pd.DataFrame) -> pd.DataFrame:
    c_counts = df_combi.groupby(['vclass', 'donor'])['gene'].nunique().unstack().fillna(0).astype(int)
    p_counts = df_ppcg.groupby(['vclass', 'donor'])['gene'].nunique().unstack().fillna(0).astype(int)
    selected = set()

    c_donors = sorted(list(df_combi['donor'].unique()))
    p_donors = sorted(list(df_ppcg['donor'].unique()))
    for c_donor in c_donors:
        diffs = []
        for p_donor in p_donors:
            if p_donor in selected:
                continue
            c = c_counts[c_donor]
            p = p_counts[p_donor]
            d = c-p 
            d = d.abs()
            diff = d.sum()
            diffs.append((p_donor, diff))
        diffs = sorted(diffs, key=lambda x: x[1])
        sel = diffs[0][0]
        selected.add(sel)
    
    df = df_ppcg[df_ppcg['donor'].isin(selected)].copy()
    return df

def load_cmdline_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Runs downstream mutation analysis.')
    
    # run ID
    parser.add_argument('--run-id', type=str, default='default', help='ID for this analysis.')
    parser.add_argument('--outdir', type=str, default='./results', help='Path to output directory.')

    # files
    parser.add_argument('--muts', type=str, default=INFILE_MUTATIONS, help='Path to assigned mutations file.')
    parser.add_argument('--sheet', type=str, default=INFILE_SHEET, help='Path to samplesheet (metadata).')
    parser.add_argument('--dpclust-dir', type=str, default=INDIR_DPCLUST, help='Path to DPClust folder')
    parser.add_argument('--conipher-dir', type=str, default=INDIR_CONIPHER, help='Path to conipher folder')
    
    # donor selection
    parser.add_argument('--ppcg-all', action='store_true', help="Use all PPCG donors rather than the 65 donor control group")
    parser.add_argument('--combi-all', action='store_true', help="Use all COMBI donors rather than the 65 with good trees")
    parser.add_argument('--hypermutator-ngenes', type=int, default=500, help="Number of mutated genes defining a hypermutator")

    # mutation selection
    parser.add_argument('--filter-secondary-leaves', action='store_true', help="Remove private secondary mutations in COMBI donors")
    parser.add_argument('--filter-primary-leaves', action='store_true', help="Remove private primary mutations in COMBI donors")
    
    # gene selection
    parser.add_argument('--top-genes', type=int, default=0, help="Only keep the <--top-genes> most commonly mutated genes")
    parser.add_argument('--min-cohort-prop', type=float, default=0.0, help="Only keep genes mutated above this proportion across both cohorts")
    parser.add_argument('--min-combi-prop', type=float, default=0.0, help="Only keep genes mutated above this proportion in combimets cohort")


    args = parser.parse_args()
    return args


if __name__ == '__main__':
    main()
