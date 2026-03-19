
import pandas as pd
import numpy as np
import gseapy as gp
import warnings

from scipy.stats import mannwhitneyu
from scipy.stats import fisher_exact, chi2
from statsmodels.formula.api import logit
from statsmodels.stats.multitest import multipletests

from pathlib import Path


# =============================================================================
# Per-gene tests (unchanged from previous script)
# =============================================================================

def run_fisher_gene_association(df: pd.DataFrame,
                                burden_col: str = "burden",
                                label_col: str = "metastatic",
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


def run_logit_gene_association(df: pd.DataFrame,
                                    burden_col: str = "burden",
                                    label_col: str = "metastatic",
                                    fdr_threshold: float = 0.05) -> pd.DataFrame:
    """
    Test each gene for association with a binary outcome using logistic regression,
    adjusting for per-patient mutational burden (total mutations excluding the
    gene being tested).

    Model per gene:
        metastatic ~ gene_mutated + burden_excluding_gene

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
        Columns: gene, n_mutated_meta, n_mutated_nonmeta, log_odds_ratio,
                 odds_ratio, p_value, q_value, significant
    """
    gene_cols = [c for c in df.columns if c != label_col and c != burden_col]
    total_muts = df[burden_col]
    meta       = df[label_col] == 1
    non_meta   = df[label_col] == 0

    results = []
    for gene in gene_cols:
        burden = total_muts - df[gene]
        burden_std = burden.std()
        burden_scaled = (burden - burden.mean()) / burden_std if burden_std > 0 else burden - burden.mean()
        mut = df[gene]
        a = (meta     & (mut == 1)).sum()
        c = (non_meta & (mut == 1)).sum()

        if mut.nunique() < 2:
            results.append({
                "gene": gene, "n_mutated_meta": int(a), "n_mutated_nonmeta": int(c),
                "log_odds_ratio": np.nan, "odds_ratio": np.nan, "p_value": np.nan,
            })
            continue

        fit_df = pd.DataFrame({
            label_col:  df[label_col].values,
            "gene_mut": mut.values,
            "burden":   burden_scaled.values,
        })
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                model   = logit(f"{label_col} ~ gene_mut + burden", data=fit_df).fit(
                    disp=False, maxiter=200)
            log_or  = model.params["gene_mut"]
            p_value = model.pvalues["gene_mut"]
        except Exception:
            log_or, p_value = np.nan, np.nan

        results.append({
            "gene":              gene,
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


# =============================================================================
# Fisher's combined p-value gene-set analysis
# =============================================================================

def run_geneset_combined(
        gene_results: pd.DataFrame,
        gene_sets: dict[str, list[str]],
        p_col: str = "p_value",
        gene_col: str = "gene",
        min_overlap: int = 3,
        n_top_genes: int = 5,
        fdr_threshold: float = 0.05,
) -> pd.DataFrame:
    """
    Gene-set association test using Fisher's combined p-value method.

    For each gene set, the per-gene p-values of its member genes are combined
    using Fisher's method:

        X = -2 * sum(ln(p_i))   ~   chi-squared(2k)

    where k is the number of genes with a valid p-value in the set.

    Parameters
    ----------
    gene_results : pd.DataFrame
        Output of run_fisher_gene_association() or run_burden_adjusted_association().
        Must contain gene_col and p_col.
    gene_sets : dict[str, list[str]]
        Gene set name -> list of gene names (e.g. from load_gmt()).
        Gene names should be upper-case to match the dataframe.
    p_col : str
        Column in gene_results containing per-gene p-values.
    gene_col : str
        Column in gene_results containing gene names.
    min_overlap : int
        Minimum number of genes that must overlap between the set and the data
        for the set to be tested. Sets below this threshold are still reported
        but with combined_p = NaN.
    n_top_genes : int
        Number of top-contributing genes (lowest p-value) to list per set.
    fdr_threshold : float
        FDR threshold for gene-set-level BH correction.

    Returns
    -------
    pd.DataFrame
        One row per gene set, sorted by combined_p, with columns:
        gene_set, set_size, n_in_data, n_tested, pct_coverage,
        combined_statistic, combined_p, q_value, significant,
        top_genes (comma-separated), missing_genes (count)
    """
    # Normalise gene names to upper-case for matching
    gene_results = gene_results.copy()
    gene_results[gene_col] = gene_results[gene_col].str.upper()

    # Index p-values by gene name for fast lookup
    p_series = (
        gene_results
        .dropna(subset=[p_col])
        .set_index(gene_col)[p_col]
    )

    results = []

    for set_name, set_genes in gene_sets.items():
        set_genes_upper = [g.upper() for g in set_genes]
        set_size        = len(set_genes_upper)

        # Genes in set that appear in the data (with a valid p-value)
        overlap_genes   = [g for g in set_genes_upper if g in p_series.index]
        n_in_data       = len(overlap_genes)
        missing_count   = set_size - n_in_data
        pct_coverage    = round(100 * n_in_data / set_size, 1) if set_size > 0 else 0.0

        # Collect p-values for overlapping genes
        subset_p = p_series.loc[overlap_genes]
        n_tested = len(subset_p)

        if n_tested < min_overlap:
            # Still report but mark as not tested
            results.append({
                "gene_set":           set_name,
                "set_size":           set_size,
                "n_in_data":          n_in_data,
                "n_tested":           n_tested,
                "pct_coverage":       pct_coverage,
                "missing_genes":      missing_count,
                "combined_statistic": np.nan,
                "combined_p":         np.nan,
                "top_genes":          "",
            })
            continue

        # Fisher's combination statistic: -2 * sum(ln(p))
        # Clip p-values away from 0 to avoid ln(0) = -inf
        p_clipped  = np.clip(subset_p.values, 1e-300, 1.0)
        statistic  = -2.0 * np.sum(np.log(p_clipped))
        combined_p = chi2.sf(statistic, df=2 * n_tested)

        # Top contributing genes = lowest individual p-values
        top = subset_p.nsmallest(n_top_genes)
        top_genes_str = ", ".join(
            f"{g}({p:.2e})" for g, p in top.items()
        )

        results.append({
            "gene_set":           set_name,
            "set_size":           set_size,
            "n_in_data":          n_in_data,
            "n_tested":           n_tested,
            "pct_coverage":       pct_coverage,
            "missing_genes":      missing_count,
            "combined_statistic": round(statistic, 4),
            "combined_p":         combined_p,
            "top_genes":          top_genes_str,
        })

    results_df = pd.DataFrame(results)

    # BH correction over testable sets only
    testable = results_df["combined_p"].notna()
    reject   = np.zeros(len(results_df), dtype=bool)
    q_vals   = np.full(len(results_df), np.nan)

    if testable.sum() > 0:
        rej_t, q_t, _, _ = multipletests(
            results_df.loc[testable, "combined_p"],
            method="fdr_bh",
            alpha=fdr_threshold,
        )
        reject[testable] = rej_t
        q_vals[testable] = q_t

    results_df["q_value"]     = q_vals
    results_df["significant"] = reject

    results_df = results_df.sort_values("combined_p").reset_index(drop=True)

    return results_df


# =============================================================================
# Summary helper
# =============================================================================

def summarise_genesets(results_df: pd.DataFrame,
                       fdr_threshold: float = 0.05,
                       n: int = 20) -> None:
    """Print a summary of the gene-set association results."""
    testable = results_df["combined_p"].notna()
    sig      = results_df[results_df["significant"]]

    print(f"Gene sets loaded:    {len(results_df)}")
    print(f"Gene sets tested:    {testable.sum()}")
    print(f"Significant (q<{fdr_threshold}): {len(sig)}")

    cols = ["gene_set", "set_size", "n_in_data", "pct_coverage",
            "combined_p", "q_value", "top_genes"]

    print(f"\nTop {n} gene sets by combined p-value:")
    top = results_df[testable].head(n)
    print(top[cols].to_string(index=False))








# =============================================================================
# ssGSEA scoring
# =============================================================================
 
def run_ssgsea(expr_matrix: pd.DataFrame,
               gene_sets: dict[str, list[str]] | str,
               min_size: int = 5,
               max_size: int = 2000,
               sample_norm: str = "rank",
               threads: int = 4) -> pd.DataFrame:
    """
    Compute per-patient ssGSEA enrichment scores for each gene set.
 
    ssGSEA (Barbie et al. 2009) estimates, for each sample independently,
    how strongly the genes in a set are enriched at the top of the
    genome-wide mutation rank. Scores are comparable within a sample;
    the resulting matrix (gene sets × patients) can be tested across
    groups with standard statistics.
 
    Parameters
    ----------
    expr_matrix : pd.DataFrame
        Genes × patients matrix (output of prepare_expression_matrix).
    gene_sets : dict or str
        Either a dict {set_name: [gene, ...]} or a path to a .gmt file.
    min_size : int
        Minimum number of genes in a set that must overlap with the data.
        Sets below this threshold are silently skipped by GSEApy.
    max_size : int
        Maximum gene set size to consider.
    sample_norm : str
        Normalisation method for GSEApy ssGSEA. 'rank' (default) converts
        expression values to ranks before scoring — appropriate for binary
        mutation data. Use 'custom' to pass raw values unchanged.
    threads : int
        Number of parallel threads for GSEApy.
 
    Returns
    -------
    pd.DataFrame
        Gene sets × patients enrichment score matrix.
        Index = gene set names, columns = patient IDs.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = gp.ssgsea(
            data              = expr_matrix,
            gene_sets         = gene_sets,
            sample_norm_method= sample_norm,
            min_size          = min_size,
            max_size          = max_size,
            threads           = threads,
            outdir            = None,
            no_plot           = True,
            verbose           = False,
        )
 
    # GSEApy returns a result object; extract the score matrix
    scores = result.res2d.pivot(index="Term", columns="Name", values="NES")
    scores.index.name   = "gene_set"
    scores.columns.name = None
    return scores
 
 
# =============================================================================
# Differential ssGSEA score test
# =============================================================================
 
def test_ssgsea_scores(scores: pd.DataFrame,
                       labels: pd.Series,
                       fdr_threshold: float = 0.05) -> pd.DataFrame:
    """
    Test each gene set for differential enrichment between metastatic and
    non-metastatic patients using the Mann-Whitney U test (Wilcoxon rank-sum).
 
    Mann-Whitney is preferred over a t-test here because ssGSEA scores are
    not normally distributed, and the test is robust to the heavy-tailed,
    sometimes bimodal distributions that arise with binary input data.
 
    The effect size reported is the rank-biserial correlation (r), a
    standardised, interpretable measure:
        r = 1   → all metastatic scores exceed all non-metastatic
        r = 0   → groups fully overlap
        r = -1  → all non-metastatic scores exceed all metastatic
 
    Parameters
    ----------
    scores : pd.DataFrame
        Gene sets × patients enrichment score matrix (from run_ssgsea).
    labels : pd.Series
        Binary outcome series indexed by patient ID (1 = metastatic).
    fdr_threshold : float
        FDR threshold for BH correction.
 
    Returns
    -------
    pd.DataFrame
        One row per gene set, sorted by p-value, with columns:
        gene_set, mean_score_meta, mean_score_nonmeta, effect_size_r,
        p_value, q_value, significant
    """
    meta_ids     = labels[labels == 1].index
    non_meta_ids = labels[labels == 0].index
 
    # Align: keep only patients present in the score matrix
    meta_ids     = [p for p in meta_ids     if p in scores.columns]
    non_meta_ids = [p for p in non_meta_ids if p in scores.columns]
 
    n_meta     = len(meta_ids)
    n_non_meta = len(non_meta_ids)
 
    results = []
    for gene_set, row in scores.iterrows():
        s_meta     = row[meta_ids].dropna().astype(float).values
        s_non_meta = row[non_meta_ids].dropna().astype(float).values
 
        if len(s_meta) < 2 or len(s_non_meta) < 2:
            results.append({
                "gene_set":          gene_set,
                "mean_score_meta":   np.nan,
                "mean_score_nonmeta":np.nan,
                "effect_size_r":     np.nan,
                "p_value":           np.nan,
            })
            continue
        
        stat, p = mannwhitneyu(s_meta, s_non_meta, alternative="two-sided")
 
        # Rank-biserial correlation as effect size
        # r = 2U / (n1 * n2) - 1  (ranges from -1 to +1)
        r = (2 * stat) / (len(s_meta) * len(s_non_meta)) - 1
 
        results.append({
            "gene_set":           gene_set,
            "mean_score_meta":    round(float(np.mean(s_meta)),     4),
            "mean_score_nonmeta": round(float(np.mean(s_non_meta)), 4),
            "effect_size_r":      round(float(r),                   4),
            "p_value":            p,
        })
 
    results_df = pd.DataFrame(results)
 
    # BH correction
    valid  = results_df["p_value"].notna()
    reject = np.zeros(len(results_df), dtype=bool)
    q_vals = np.full(len(results_df), np.nan)
 
    if valid.sum() > 0:
        rej_v, q_v, _, _ = multipletests(
            results_df.loc[valid, "p_value"], method="fdr_bh", alpha=fdr_threshold
        )
        reject[valid] = rej_v
        q_vals[valid] = q_v
 
    results_df["q_value"]     = q_vals
    results_df["significant"] = reject
 
    results_df = (
        results_df
        .sort_values("p_value")
        .reset_index(drop=True)
    )
 
    return results_df
 


#####################
### GENESET LOGIT ###
#####################

def run_logit_gset_association(
    df: pd.DataFrame,
    burden_LUT: dict[str, int],
    label_col: str = "metastatic",
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
    gset_cols  = [c for c in df.columns if c != label_col]
    total_muts = pd.Series(burden_LUT)
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

