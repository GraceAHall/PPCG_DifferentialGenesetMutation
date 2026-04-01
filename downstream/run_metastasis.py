
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
from filtering import filter_hypermutators
from selection import select_genesets

from reporting import (
    summarise_basic_info,
    summarise_vclasses,
    summarise_annotations,
    summarise_donors,
    plot_distribution
)

from formatting import generate_geneset_matrix


# from formatting import generate_gene_matrix
# from formatting import generate_gene_matrix_4genesets
# from formatting import prepare_gsea_matrix
# from formatting import generate_geneset_matrix

# from stats import run_fisher_gene_association
# from stats import run_logit_gene_association
# from stats import run_logit_gset_association
# from stats import run_geneset_combined
# from stats import run_ssgsea
# from stats import test_ssgsea_scores

import warnings
warnings.filterwarnings('ignore')

pd.options.display.float_format = "{:,.3f}".format
pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 200)
pd.set_option('display.width', 180)

# INFILE_MUTATIONS = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/manual/mutations.assigned.160326.tsv'
# INFILE_SIZES = '/home/grace/work/PPCG_DifferentialGenesetMutation/outputs/angel_mutations_all_donors/variant_processing/sizes.tsv'
# INFILE_SHEET = '/home/grace/work/PPCG_DifferentialGenesetMutation/samplesheet.angel.alldonors.tsv'
# INFILE_HGNC = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/supporting/hgnc/hgnc_complete_set.txt'
# INFILE_GENESETS = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/genesets/c2.cp.reactome.v2026.1.Hs.symbols.gmt'
# INFILE_GENESETS = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/genesets/h.all.v2026.1.Hs.symbols.gmt'
# INFILE_GENESETS = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/genesets/c5.go.bp.v2026.1.Hs.symbols.gmt'
# INFILE_GENESETS = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/genesets/c2.cp.kegg_medicus.v2026.1.Hs.symbols.gmt'
# INDIR_DPCLUST = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/angel/dpclust_ccf'
# INDIR_CONIPHER = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/angel/conipher_trees'


def main():
    args = load_cmdline_args()
    outdir = f"{args.outdir}/{args.run_id}"
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    report_settings(args)

    # mutations
    muts = load_mutations(args.muts, args.sheet)

    # filtering
    muts = filter_mutations(muts, args)

    # reporting
    report_mutations(muts, args)

    # genesets
    gset_LUT = load_gmt(args.gmt)
    gset_LUT, sel = select_genesets(gset_LUT, muts)
    
    # enrichment
    fisher_association(muts, gset_LUT, args)
    sys.exit()
    logistic_association(muts, gset_LUT, args)
    background_enrichment(muts, gset_LUT, args)


#########################
### STATISTICAL TESTS ###
#########################

def background_enrichment(muts: pd.DataFrame, gset_LUT: dict, args: argparse.Namespace) -> None:
    raise NotImplementedError
    # donors = set(muts[muts['cohort']=='COMBI']['donor'].unique())
    # seq_muts = muts[(muts['cohort']=='COMBI') & (muts['vclass'].isin(['SNV', 'INDEL']))].copy()
    # struct_muts = muts[(muts['cohort']=='COMBI') & (muts['vclass']=='SV')].copy()
    # seqmat = generate_gene_matrix_4genesets(seq_muts, donors, args)
    # structmat = generate_gene_matrix_4genesets(struct_muts, donors, args)
    # print(seqmat.shape)
    # print(structmat.shape)
    # do_geneset_enrichment_bootstrapping(gset_LUT, seqmat, structmat, args)

def logistic_association(muts: pd.DataFrame, gset_LUT: dict, args: argparse.Namespace) -> None:
    all_donors = sorted(list(muts['donor'].unique()))
    mat = generate_geneset_matrix(gset_LUT, muts, all_donors)
    raise NotImplementedError
    # mat = generate_geneset_matrix(gset_LUT, muts, args)
    # burden_LUT = muts.groupby('donor')['gene'].nunique().to_dict()
    # res_logit = run_logit_gset_association(mat, burden_LUT=burden_LUT)
    # res_logit = res_logit.sort_values('p_value')
    # print('\nLOGIT')
    # print(res_logit.shape)
    # print(res_logit.head(10))
    # outfile_logit = f"{args.outdir}/{args.run_id}/logit.tsv"
    # res_logit.to_csv(outfile_logit, sep='\t', index=False, float_format='%.5f')

def fisher_association(muts: pd.DataFrame, gset_LUT: dict, args: argparse.Namespace) -> None:
    res_fish = _do_fisher_association(gset_LUT, muts)
    res_fish = res_fish.sort_values('p_value')
    print('\nFISHER')
    print(res_fish.shape)
    print(res_fish.head(10))
    outfile_fish = f"{args.outdir}/{args.run_id}/fisher.tsv"
    res_fish.to_csv(outfile_fish, sep='\t', index=False, float_format='%.5f')

def _do_fisher_association(gset_LUT: dict[str, list[str]], muts: pd.DataFrame) -> pd.DataFrame:
    
    # Define a function to apply the test to each row
    def apply_fisher_test(row):
        table = [[row['a'], row['b']], [row['c'], row['d']]]
        odds_ratio, p_value = fisher_exact(table, alternative='two-sided')
        return pd.Series({'odds_ratio': odds_ratio, 'p_value': p_value})
    
    data = []
    n_donors_combi = muts[muts['cohort']=='COMBI']['donor'].nunique()
    n_donors_ppcg = muts[muts['cohort']=='PPCG']['donor'].nunique()
    for gset, genes in gset_LUT.items():
        df = muts[muts['gene'].isin(set(genes))]
        n_genes_combi = df[df['cohort']=='COMBI']['gene'].nunique()
        n_genes_ppcg = df[df['cohort']=='PPCG']['gene'].nunique()
        a = df[df['cohort']=='COMBI']['donor'].nunique()
        c = df[df['cohort']=='PPCG']['donor'].nunique()
        data.append((gset, len(genes), n_genes_combi, n_genes_ppcg, a, n_donors_combi-a, c, n_donors_ppcg-c))
    df = pd.DataFrame.from_records(data, columns=['geneset', 'total_genes', 'mut_genes_combi', 'mut_genes_ppcg', 'a', 'b', 'c', 'd'])

    # Apply the function across rows (axis=1) and add results as new columns
    res = df.apply(apply_fisher_test, axis=1)
    df = pd.concat([df, res], axis=1)
    df = df.sort_values('p_value')
    df = df.set_index('geneset')

    reject, q_values, _, _ = multipletests(
        df["p_value"], method="fdr_bh", alpha=0.05
    )
    df["q_value"]     = q_values
    df["significant"] = reject
    df = df.reset_index()

    df = df[[x for x in df.columns if x!='geneset']+['geneset']].copy()
    return df 

#################
### REPORTING ###
#################

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
        summarise_annotations(df)
        summarise_donors(df)
    sys.stdout = original_stdout

    print(f'writing plot to {outfile_plot}')
    plot_distribution(df, filepath=outfile_plot)


################################################################
### LOADING MUTATIONS, FILTERING MUTATIONS, CREATING COHORTS ###
################################################################

def filter_mutations(df: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    print('filtering donors/mutations')
    df = df[df['annotation']!='LOH'].copy()
    df_ppcg = df[df['cohort']=='PPCG'].copy()
    df_combi = df[df['cohort']=='COMBI'].copy()
    df_combi = filter_seedtraj_donors(df_combi, args)
    df_combi = filter_seedtraj_mutations(df_combi, args)
    df_combi = filter_hypermutators(df_combi, max_mutated_genes=args.hypermutator_ngenes)
    df_ppcg = filter_hypermutators(df_ppcg, max_mutated_genes=args.hypermutator_ngenes)
    df_ppcg = pick_matched_cohort(df_ppcg, df_combi)

    df = pd.concat([df_ppcg, df_combi], ignore_index=True)
    return df 

def filter_seedtraj_donors(table: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    df = table.copy()
    meta = pd.read_csv(args.clone_meta, sep='\t', header=0)
    target_donors = set(meta['patient'].unique())
    df = df[df['donor'].isin(target_donors)].copy()
    return df 
    
def filter_seedtraj_mutations(table: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    df = table.copy()

    donor2clones = dict()
    meta = pd.read_csv(args.clone_meta, sep='\t', header=0)
    meta['clone'] = meta['clone'].astype(str)
    for donor, donor_df in meta.groupby('patient'):
        # get clones considered on seeding trajectory
        path_clones = set()
        for traj in donor_df['ancestors'].values:
            nodelist = traj.split('-')
            nodelist = nodelist[:-1]
            assert len(nodelist) >= 1
            path_clones.update(nodelist)
        donor2clones[donor] = path_clones

    df['seeding_trajectory'] = False 
    for donor in df['donor'].unique():
        mask = (df['donor']==donor) & (df['clone'].isin(donor2clones[donor]))
        df.loc[mask, 'seeding_trajectory'] = True 
    df = df[df['seeding_trajectory']==True].copy()
    df = df.drop('seeding_trajectory', axis=1)
    return df.copy()

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
        # print('\nSELECTED')
        # print()
        # print(c_counts[c_donor])
        # print()
        # print(p_counts[sel])
        selected.add(sel)
    
    df = df_ppcg[df_ppcg['donor'].isin(selected)].copy()
    return df




def do_geneset_enrichment_bootstrapping(
    gset_LUT: dict[str, list[str]], 
    seqmat: pd.DataFrame, 
    structmat: pd.DataFrame, 
    args: argparse.Namespace) -> pd.DataFrame:
    
    print('\nBOOTSTRAPPING')
    BUDGET = 1000

    # gene sizes for seqvars and structvars 
    sizes = pd.read_csv(args.sizes, sep='\t', header=0)
    seq_sizes_LUT = sizes.set_index('gene')['cum_exon_len'].astype(float).to_dict()
    struct_sizes_LUT = sizes.set_index('gene')['span'].astype(float).to_dict()
    
    # prep    
    seqmat_obs = seqmat.T
    structmat_obs = structmat.T
    
    # observed
    mat_obs = seqmat_obs + structmat_obs
    mat_obs = mat_obs.clip(upper=1)
    res_obs = count_mutated_donors(mat_obs, gset_LUT)
    
    # random permutations
    perms = {gset: [] for gset in list(gset_LUT.keys())}
    for i in range(BUDGET):
        if i % 10 == 0:
            print(i)
        
        seqmat_p = permute_mutation_matrix_weighted(seqmat_obs, seq_sizes_LUT)
        structmat_p = permute_mutation_matrix_weighted(structmat_obs, struct_sizes_LUT)
        mat_p = seqmat_p + structmat_p
        mat_p = mat_p.clip(upper=1)
        res_p = count_mutated_donors(mat_p, gset_LUT)
        for gset, count in res_p.to_dict().items():
            perms[gset].append(count)
    
    # results
    sframe = pd.DataFrame(index=res_obs.index.to_list())
    sframe['observed'] = res_obs
    sframe['background'] = perms
    sframe = test_enrichment_empirical(sframe)
    sframe['background'] = sframe['background'].apply(lambda x: sum(x)/len(x))

    # final formatting
    res = sframe.copy()
    print(res['significant'].value_counts())
    res = res.reset_index()
    res = res.rename(columns={'index': 'geneset'})
    res = res[[x for x in res.columns if x!='geneset']+['geneset']].copy()
    outfile_enrich = f"{args.outdir}/{args.run_id}/enrichment.tsv"
    res.to_csv(outfile_enrich, sep='\t', index=False, float_format='%.5f')
    return res 




def test_enrichment_empirical(df: pd.DataFrame,
                               observed_col: str = "observed",
                               background_col: str = "background",
                               fdr_threshold: float = 0.05) -> pd.DataFrame:
    """
    Test whether each gene set is mutated significantly above its bootstrap
    background using an empirical (permutation) p-value.
 
    For each gene set the empirical p-value is:
 
        p = (# replicates >= observed + 1) / (N + 1)
 
    The +1 pseudocount in numerator and denominator is the standard correction
    that avoids p = 0 when the observed value exceeds every replicate, which
    would be an overconfident claim given finite bootstrap sampling.
 
    Because this is a one-sided test (enrichment above background), p-values
    are not doubled. If you also want to test for depletion, flip the
    inequality to <= and run separately.
 
    Multiple testing correction is applied across gene sets using
    Benjamini-Hochberg FDR.
 
    Parameters
    ----------
    df : pd.DataFrame
        Rows = gene sets. Must contain:
            observed_col   : int, observed mutated donor count.
            background_col : list of int/float, one value per bootstrap
                             replicate.
    observed_col : str
        Column name for the observed counts.
    background_col : str
        Column name containing the bootstrap distributions (lists).
    fdr_threshold : float
        FDR significance threshold for BH correction.
 
    Returns
    -------
    pd.DataFrame
        Original dataframe with four new columns appended:
            n_replicates  : number of bootstrap replicates used
            p_value       : empirical one-sided p-value
            q_value       : BH-corrected FDR
            significant   : bool, True if q_value < fdr_threshold
    """
    result = df.copy()
 
    p_values     = []
    n_replicates = []
 
    for _, row in df.iterrows():
        observed   = row[observed_col]
        background = np.asarray(row[background_col])
        n          = len(background)
 
        p = (np.sum(background >= observed) + 1) / (n + 1)
 
        p_values.append(p)
        n_replicates.append(n)
 
    result["n_replicates"] = n_replicates
    result["p_value"]      = p_values
 
    reject, q_values, _, _ = multipletests(p_values, method="fdr_bh",
                                           alpha=fdr_threshold)
    result["q_value"]     = q_values
    result["significant"] = reject
 
    return result.sort_values("p_value")
 
 
    # for gset in pframe.columns:
    #     obs = int(sframe.loc[gset, 'observed'])
    #     n_perms_gteq = (pframe.T[gset]>=obs).sum()
    #     sframe.loc[gset, 'p-value']

    # pframe['p-value'] = pframe.apply(lambda x: x>=, axis=1)



    # permcols = [x for x in pframe.cols if x.startswith('perm')]
    

    # print(sframe.iloc[:, 1:].head())
    # sframe['odds_ratio'] = sframe['observed'] - sframe['background']
    # print()
    # print(sframe.head(10))
    # print()
    # sframe = sframe.sort_values('diff', ascending=False)
    # sframe = sframe
    # print(sframe.head(10))


def count_mutated_donors(mutation_matrix: pd.DataFrame,
                         gene_sets: dict[str, list[str]]) -> pd.Series:
    """
    For each gene set, count the number of donors (patients) who have at
    least one mutation in any gene belonging to that set.
 
    Parameters
    ----------
    mutation_matrix : pd.DataFrame
        Binary matrix with genes as rows and patients as columns.
        Values must be 0 or 1.
    gene_sets : dict[str, list[str]]
        Mapping of gene set name -> list of gene names.
        Gene names that are absent from the mutation matrix are silently
        ignored (only the overlap is used).
 
    Returns
    -------
    pd.Series
        Integer series indexed by gene set name, giving the number of
        donors with at least one mutation in that set.
        Gene sets with no overlapping genes return 0.
 
    Examples
    --------
    >>> mm = pd.DataFrame(
    ...     [[1, 0], [0, 1], [0, 0]],
    ...     index=["BRCA1", "TP53", "EGFR"],
    ...     columns=["P1", "P2"],
    ... )
    >>> gs = {"set_A": ["BRCA1", "EGFR"], "set_B": ["TP53", "EGFR"]}
    >>> count_mutated_donors(mm, gs)
    set_A    1
    set_B    1
    dtype: int64
    """
    counts = {}
    for set_name, genes in gene_sets.items():
        overlap = [g for g in genes if g in mutation_matrix.index]
        if overlap:
            # any() across rows: True if any gene in set is mutated per patient
            mutated = mutation_matrix.loc[overlap].any(axis=0)
            counts[set_name] = int(mutated.sum())
        else:
            counts[set_name] = 0
 
    return pd.Series(counts, dtype=int)
 
 
def permute_mutation_matrix_weighted(mutation_matrix: pd.DataFrame,
                             gene_sizes: dict[str, float] | None = None,
                             rng: np.random.Generator | None = None) -> pd.DataFrame:
    """
    Return a randomly permuted copy of a binary mutation matrix in which
    each patient's total mutation count (column sum) is exactly preserved.

    For each patient (column), the mutations are redistributed across genes
    by weighted sampling without replacement. Row sums (per-gene mutation
    frequencies) are NOT preserved — genes are sampled freely.

    Parameters
    ----------
    mutation_matrix : pd.DataFrame
        Binary matrix with genes as rows and patients as columns.
        Values must be 0 or 1.
    gene_sizes : dict[str, float], optional
        Mapping of gene name -> size (e.g. coding length in base pairs).
        Sampling probability for each gene is proportional to its size, so
        a gene of size 1000 is sampled 10x more often than one of size 100.
        Genes absent from this dictionary are assigned the median size of
        all provided values. If None, all genes are sampled with equal
        probability (uniform permutation).
    rng : np.random.Generator, optional
        A numpy random Generator for reproducibility, e.g.
        np.random.default_rng(seed=42). If None, a fresh Generator is
        created using entropy from the OS.

    Returns
    -------
    pd.DataFrame
        Permuted binary matrix with the same shape, index, and columns as
        the input. Each column sum equals the corresponding column sum in
        the original matrix.

    Examples
    --------
    >>> rng = np.random.default_rng(42)
    >>> df = pd.DataFrame([[1,0],[0,1],[1,1],[0,0]], columns=["P1","P2"])
    >>> permuted = permute_mutation_matrix(df, rng=rng)
    >>> (permuted.sum(axis=0) == df.sum(axis=0)).all()
    True
    """
    if rng is None:
        rng = np.random.default_rng()

    n_genes  = len(mutation_matrix)
    col_sums = mutation_matrix.sum(axis=0).astype(int)
    permuted = np.zeros((n_genes, len(mutation_matrix.columns)), dtype=np.int8)

    # Build normalised sampling weights aligned to the matrix row order
    if gene_sizes is not None:
        provided_sizes = np.array(list(gene_sizes.values()), dtype=float)
        median_size    = float(np.median(provided_sizes))

        raw_weights = np.array(
            [gene_sizes.get(g, median_size) for g in mutation_matrix.index],
            dtype=float,
        )
        # Normalise to a probability distribution
        sampling_probs = raw_weights / raw_weights.sum()
    else:
        sampling_probs = None   # rng.choice uses uniform by default

    for j, n_mutations in enumerate(col_sums):
        if n_mutations > 0:
            mutated_rows = rng.choice(n_genes, size=n_mutations,
                                      replace=False, p=sampling_probs)
            permuted[mutated_rows, j] = 1

    return pd.DataFrame(
        permuted,
        index   = mutation_matrix.index,
        columns = mutation_matrix.columns,
    )



def permute_mutation_matrix_uniform(mutation_matrix: pd.DataFrame,
                             rng: np.random.Generator | None = None) -> pd.DataFrame:
    """
    Return a randomly permuted copy of a binary mutation matrix in which
    each patient's total mutation count (column sum) is exactly preserved.
 
    For each patient (column), the mutations are redistributed uniformly at
    random across all genes, without replacement. Row sums (per-gene mutation
    frequencies) are NOT preserved — genes are sampled freely.
 
    Parameters
    ----------
    mutation_matrix : pd.DataFrame
        Binary matrix with genes as rows and patients as columns.
        Values must be 0 or 1.
    rng : np.random.Generator, optional
        A numpy random Generator for reproducibility, e.g.
        np.random.default_rng(seed=42). If None, a fresh Generator is
        created using entropy from the OS.
 
    Returns
    -------
    pd.DataFrame
        Permuted binary matrix with the same shape, index, and columns as
        the input. Each column sum equals the corresponding column sum in
        the original matrix.
 
    Examples
    --------
    >>> rng = np.random.default_rng(42)
    >>> df = pd.DataFrame([[1,0],[0,1],[1,1],[0,0]], columns=["P1","P2"])
    >>> permuted = permute_mutation_matrix(df, rng=rng)
    >>> (permuted.sum(axis=0) == df.sum(axis=0)).all()
    True
    """
    if rng is None:
        rng = np.random.default_rng()
 
    n_genes     = len(mutation_matrix)
    col_sums    = mutation_matrix.sum(axis=0).astype(int)
    permuted    = np.zeros((n_genes, len(mutation_matrix.columns)), dtype=np.int8)
 
    for j, n_mutations in enumerate(col_sums):
        if n_mutations > 0:
            mutated_rows = rng.choice(n_genes, size=n_mutations, replace=False)
            permuted[mutated_rows, j] = 1
 
    return pd.DataFrame(
        permuted,
        index   = mutation_matrix.index,
        columns = mutation_matrix.columns,
    )
 

def do_geneset_enrichment_gsea(muts: pd.DataFrame, mat: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    print("\nPreparing expression matrix ...")
    expr_matrix = prepare_gsea_matrix(mat)
    print()
    print(expr_matrix.iloc[:5, :5])
    print()
    print(f"  Matrix shape: {expr_matrix.shape[0]} genes × {expr_matrix.shape[1]} patients\n")

    print("Running ssGSEA ...")
    gset_LUT = load_gmt(args.gmt)
    scores = run_ssgsea(
        expr_matrix = expr_matrix,
        gene_sets   = gset_LUT,
        min_size    = 5,
        threads     = 4,
    )
    print()
    print(scores.iloc[:5, :5])
    print()
    print(f"  Score matrix: {scores.shape[0]} gene sets × {scores.shape[1]} patients\n")

    print("Testing ssGSEA scores (Mann-Whitney U) ...")
    labels  = mat["metastatic"]
    results = test_ssgsea_scores(scores, labels, fdr_threshold=0.05)

    counts_combi = {}
    counts_ppcg = {}
    for gset, genes in gset_LUT.items():
        df = muts[muts['gene'].isin(set(genes))]
        counts_combi[gset] = df[df['cohort']=='COMBI']['donor'].nunique()
        counts_ppcg[gset] = df[df['cohort']=='PPCG']['donor'].nunique()
    results = results.set_index('gene_set')
    results['COMBI'] = counts_combi
    results['PPCG'] = counts_ppcg
    n_combi = muts[muts['cohort']=='COMBI']['donor'].nunique()
    n_ppcg = muts[muts['cohort']=='PPCG']['donor'].nunique()
    results['odds_ratio'] = (results['COMBI']/n_combi) / (results['PPCG']/n_ppcg)
    results = results.reset_index()
    print()
    print(results[results['COMBI']>results['PPCG']].head(10))
    print()
    print(results['significant'].value_counts())
    # summarise(results, fdr_threshold=0.05, n=10)

    import seaborn as sns 
    import matplotlib.pyplot as plt 
    outfile_clustermap = f"{args.outdir}/{args.run_id}/clustermap.png"
    targets = results.head(50)['gene_set'].to_list()
    # donor2cohort = muts.drop_duplicates('donor').set_index('donor')['cohort'].to_dict()
    # sns.clustermap(scores.loc[targets])

    donors_combi = set(muts[muts['cohort']=='COMBI']['donor'].unique())
    col_colors = ['red' if x in donors_combi else 'blue' for x in scores.columns]
    sns.clustermap(data=scores.loc[targets].astype(float), col_colors=col_colors, figsize=(40, 20))
    plt.savefig(outfile_clustermap, dpi=300)
    plt.close()

def do_geneset_enrichment_combp(muts: pd.DataFrame, res_gene: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    print('running gene set enrichment tests')
    outfile_gsets = f"{args.outdir}/{args.run_id}/genesets.tsv"
    gset_LUT = load_gmt(args.gmt)

    counts_combi = {}
    counts_ppcg = {}
    for gset, genes in gset_LUT.items():
        df = muts[muts['gene'].isin(set(genes))]
        counts_combi[gset] = df[df['cohort']=='COMBI']['donor'].nunique()
        counts_ppcg[gset] = df[df['cohort']=='PPCG']['donor'].nunique()
    
    res = run_geneset_combined(res_gene, gene_sets=gset_LUT)
    res = res.set_index('gene_set')
    res['COMBI'] = counts_combi
    res['PPCG'] = counts_ppcg
    res = res.reset_index()

    res = res[[x for x in res.columns if x!='gene_set'] + ['gene_set']]
    res.to_csv(outfile_gsets, sep='\t', index=False, float_format='%.4f')
    res['pct_coverage'] = res['pct_coverage'].astype(int)
    res = res.sort_values('combined_p')
    print(res.shape)
    print(res.head(10))
    return res

def do_gene_enrichment(df: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    print('running gene enrichment tests')
    outfile_fish = f"{args.outdir}/{args.run_id}/gene_fisher.tsv"
    outfile_logit = f"{args.outdir}/{args.run_id}/gene_logit.tsv"
    outfile_comb = f"{args.outdir}/{args.run_id}/gene_combined.tsv"
    
    res_fish = run_fisher_gene_association(df)
    res_logit = run_logit_gene_association(df)
    res_fish['odds_ratio'] = res_fish['odds_ratio'].clip(upper=100)
    res_logit['odds_ratio'] = res_logit['odds_ratio'].clip(upper=100)
    res_logit['log_odds_ratio'] = res_logit['log_odds_ratio'].clip(upper=100)

    res_fish.to_csv(outfile_fish, sep='\t', index=False, float_format='%.5f')
    res_logit.to_csv(outfile_logit, sep='\t', index=False, float_format='%.5f')

    res_fish = res_fish.set_index('gene')
    res_logit = res_logit.set_index('gene')
    res = res_logit.copy()
    mask = (res['n_mutated_meta']==0) | (res['n_mutated_nonmeta']==0)
    res.loc[mask, 'p_value'] = res_fish.loc[mask, 'p_value']
    res.loc[mask, 'q_value'] = res_fish.loc[mask, 'q_value']
    res.loc[mask, 'significant'] = res_fish.loc[mask, 'significant']
    res = res.reset_index()
    
    res.to_csv(outfile_comb, sep='\t', index=False, float_format='%.5f')
    return res


# def filter_mutations(df: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
#     print('filtering donors/mutations')
#     df_ppcg = df[df['cohort']=='PPCG'].copy()
#     df_combi = df[df['cohort']=='COMBI'].copy()
    
#     if not args.ppcg_all:
#         df_ppcg = filter_ppcg_control(df_ppcg)
    
#     if not args.combi_all:
#         df_combi = filter_donors_without_prostate(df_combi)
#         df_combi = filter_donors_without_dpclust(df_combi, args.dpclust_dir)
#         df_combi = filter_donors_without_trees(df_combi, args.conipher_dir)

#     if args.filter_primary_leaves:
#         df_combi = filter_primary_leaves(df_combi, verbose=True)
#     if args.filter_secondary_leaves:
#         df_combi = filter_secondary_leaves(df_combi, verbose=True)

#     df = pd.concat([df_ppcg, df_combi], ignore_index=True)
#     df = filter_hypermutators(df, max_mutated_genes=args.hypermutator_ngenes)
#     return df 



def load_cmdline_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Runs downstream mutation analysis.')
    parser.add_argument('--run-id', type=str, default='default', help='ID for this analysis.')
    parser.add_argument('--outdir', type=str, default='./results', help='Path to output directory.')
    parser.add_argument('--gmt', type=str, required=True, help='Path to msigdb genesets .gmt file.')
    parser.add_argument('--muts', type=str, required=True, help='Path to assigned mutations file.')
    parser.add_argument('--sheet', type=str, required=True, help='Path to samplesheet (metadata).')
    parser.add_argument('--sizes', type=str, required=True, help='Path to gene sizes.')
    parser.add_argument('--clone-meta', type=str, required=True, help='Path to Angel metastatic clone metadata.')
    parser.add_argument('--hypermutator-ngenes', type=int, default=700, help="Number of mutated genes defining a hypermutator")

    # parser.add_argument('--gmt', type=str, default=INFILE_GENESETS, help='Path to msigdb genesets .gmt file.')
    # parser.add_argument('--muts', type=str, default=INFILE_MUTATIONS, help='Path to assigned mutations file.')
    # parser.add_argument('--sizes', type=str, default=INFILE_SIZES, help='Path to gene sizes.')
    # parser.add_argument('--hgnc', type=str, default=INFILE_HGNC, help='Path to hgnc_complete_set.txt')
    # parser.add_argument('--dpclust-dir', type=str, default=INDIR_DPCLUST, help='Path to DPClust folder')
    # parser.add_argument('--conipher-dir', type=str, default=INDIR_CONIPHER, help='Path to conipher folder')
    
    # donor selection

    # # gene selection
    # parser.add_argument('--top-genes', type=int, default=0, help="Only keep the <--top-genes> most commonly mutated genes")
    # parser.add_argument('--min-cohort-prop', type=float, default=0.0, help="Only keep genes mutated above this proportion across both cohorts")
    # parser.add_argument('--min-combi-prop', type=float, default=0.0, help="Only keep genes mutated above this proportion in combimets cohort")


    args = parser.parse_args()
    return args


if __name__ == '__main__':
    main()
