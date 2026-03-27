
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
INDIR_DPCLUST = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/angel/dpclust_ccf'
INDIR_CONIPHER = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/angel/conipher_trees'
INDIR_TIMING = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/angel/time_trees_acceleration_5'
INFILE_GENESETS = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/genesets/c5.go.bp.v2026.1.Hs.symbols.gmt'
# INFILE_GENESETS = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/genesets/c2.cp.kegg_medicus.v2026.1.Hs.symbols.gmt'
# INFILE_GENESETS = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/genesets/c2.cp.reactome.v2026.1.Hs.symbols.gmt'
# INFILE_GENESETS = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/genesets/h.all.v2026.1.Hs.symbols.gmt'

QUERIES = [
    'GOBP_RESPONSE_TO_ANTIBIOTIC'
]

def main():
    args = load_cmdline_args()
    rundir = f"{args.outdir}/{args.run_id}"
    if not os.path.exists(rundir):
        os.makedirs(rundir, exist_ok=True)

    ### mutations
    muts = load_mutations(args.muts, args.sheet)
    muts = filter_mutations(muts, args)

    gset_LUT = load_gmt(args.gmt)
    gset_LUT, sframe = select_genesets(gset_LUT, muts) 
    print(sframe[sframe['geneset'].isin(QUERIES)])

    pfields = []
    for query in QUERIES:
        genelist = gset_LUT[query]
        df = muts[muts['gene'].isin(genelist)]
        print()
        print(query)
        print()
        print(df.groupby('gene')['donor'].nunique().sort_values(ascending=False))
        print()
        print(df)

def filter_mutations(df: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    print('filtering donors/mutations')
    df = df[df['cohort']=='COMBI'].copy()
    df = df[df['vclass']!='CNA'].copy()
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
