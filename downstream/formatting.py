
import argparse
import pandas as pd 


def generate_geneset_matrix(gset_LUT: dict[str, list[str]], table: pd.DataFrame, all_donors: list[str]) -> pd.DataFrame:
    # init dataframe
    df = pd.DataFrame(index=all_donors)

    # for each geneset, add column annotating whether the donor has a mutation or not. 
    i = 0
    df = df.reset_index()
    df = df.rename(columns={'index': 'patient'})
    for gset, genelist in gset_LUT.items():
        donors = set(table[table['gene'].isin(set(genelist))]['donor'].unique())
        df[gset] = df['patient'].apply(lambda x: 1 if x in donors else 0)
        i += 1
        if i % 50 == 0:
            df = df.copy()
    df = df.set_index('patient')
    return df

def generate_gene_matrix_4genesets(table: pd.DataFrame, donors: set[str], args: argparse.Namespace) -> pd.DataFrame:
    hgnc = pd.read_csv(args.hgnc, sep='\t', header=0)
    hgnc_genes = set(hgnc[hgnc['locus_group']=='protein-coding gene']['symbol'].unique())
    mut_genes = set(table['gene'].unique())
    all_genes = sorted(list(mut_genes | hgnc_genes))
    all_donors = sorted(list(donors))

    # init dataframe
    df = pd.DataFrame(index=all_donors)

    # for each gene, add column annotating whether the donor has a mutation or not. 
    i = 0
    df = df.reset_index()
    df = df.rename(columns={'index': 'patient'})
    for gene, donors in table.groupby('gene')['donor'].agg(set).to_dict().items():
        df[gene] = df['patient'].apply(lambda x: 1 if x in donors else 0)
        i += 1
        if i % 50 == 0:
            df = df.copy()
    for gene in all_genes:
        if gene not in df.columns:
            df[gene] = 0

    df = df.set_index('patient')
    return df

def generate_gene_matrix(table: pd.DataFrame, genes: set[str]) -> pd.DataFrame:
    print('generating mutation matrix')
    all_donors = sorted(list(table['donor'].unique()))
    all_genes = sorted(list(genes))

    burden_LUT = table.groupby('donor')['gene'].nunique().to_dict()
    cohort_LUT = table.drop_duplicates('donor').set_index('donor')['cohort'].to_dict()

    # init dataframe
    df = pd.DataFrame(index=all_donors)

    # annotate membership & TMB
    df['metastatic'] = cohort_LUT
    df['metastatic'] = df['metastatic'].map({'COMBI': 1, 'PPCG': 0})
    df['metastatic'] = df['metastatic'].astype(int)
    df['burden'] = burden_LUT
    assert df['metastatic'].isna().sum() == 0
    assert df['burden'].isna().sum() == 0

    # for each gene, add column annotating whether the donor has a mutation or not. 
    i = 0
    df = df.reset_index()
    df = df.rename(columns={'index': 'patient'})
    for gene, donors in table.groupby('gene')['donor'].agg(set).to_dict().items():
        if gene in all_genes:
            df[gene] = df['patient'].apply(lambda x: 1 if x in donors else 0)
            i += 1
            if i % 50 == 0:
                df = df.copy()

    df = df.set_index('patient')
    return df


# =============================================================================
# Data formatting
# =============================================================================
 
def prepare_gsea_matrix(df: pd.DataFrame,
                               label_col: str = "metastatic") -> pd.DataFrame:
    """
    Convert the patient × gene mutation dataframe into the genes × samples
    format expected by GSEApy.
 
    GSEApy expects:
        - rows    = genes
        - columns = sample/patient IDs
        - values  = numeric (binary 0/1 mutation calls work fine for ssGSEA)
 
    Parameters
    ----------
    df : pd.DataFrame
        Patient × gene binary mutation dataframe, index = patient_id.
    label_col : str
        Name of the outcome column to drop before transposing.
 
    Returns
    -------
    pd.DataFrame
        Genes × patients matrix.
    """
    expr = df.drop(columns=[label_col]).T
    expr.index.name = "Gene"
    return expr