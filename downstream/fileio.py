
import pandas as pd 


# =============================================================================
# GMT parser
# =============================================================================

def load_gmt(gmt_path: str) -> dict[str, list[str]]:
    """
    Parse a MSigDB .gmt file into a dict mapping gene set name -> list of genes.

    GMT format (tab-delimited):
        GENE_SET_NAME  <description/url>  GENE1  GENE2  ...
    """
    gene_sets = {}
    with open(gmt_path, "r") as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            name  = parts[0]
            genes = [g.upper() for g in parts[2:] if g]
            gene_sets[name] = genes
    return gene_sets


##################
### other shyt ###
##################

def load_mutations(mutations_path: str, samplesheet_path: str, use_dpc_asmts: bool=False) -> pd.DataFrame:
    print('loading mutations')
    df = pd.read_csv(mutations_path, sep='\t', header=0)
    df['donor'] = df['sample'].apply(lambda x: x[:8])
    if use_dpc_asmts:
        hframe = load_hframe(samplesheet_path)
        df = annotate_assignment_verdict_dynamic(df, hframe)
    else:
        df = annotate_assignment_verdict_var(df)
    return df 

def load_hframe(filepath: str) -> pd.DataFrame:
    sheet = pd.read_csv(filepath, sep='\t', header=0)
    sheet['tissue'] = sheet['tissue'].replace('Recurrence', 'Metastasis')
    sheet = sheet.sort_values('sample')

    # annotating which files are present
    vclasses = ['SNV', 'INDEL', 'CNA', 'SV']
    sheet['Files'] = sheet['Missing_Files'].apply(format_files)
    for vclass in vclasses:
        sheet[vclass] = sheet['Files'].apply(lambda x: vclass in x)
    sheet = sheet.drop(['Files'], axis=1)

    data = []
    cohort_LUT = sheet.drop_duplicates('donor').set_index('donor')['cohort'].to_dict()
    for vclass in vclasses:
        for donor in sheet['donor'].unique():
            label = f"{donor}|{vclass}"
            cohort = cohort_LUT[donor]
            temp = sheet[(sheet['donor']==donor) & (sheet[vclass]==True)]
            primary = temp[temp['tissue']=='Primary']['sample'].nunique()
            secondary = temp[temp['tissue']=='Metastasis']['sample'].nunique()
            data.append((label, donor, vclass, cohort, primary, secondary))
    hframe = pd.DataFrame(data, columns=['label', 'donor', 'vclass', 'cohort', 'primary', 'secondary'])
    hframe['handling'] = hframe.apply(get_handling, axis=1)
    return hframe

def format_files(missing_text: str) -> set[str]:
    expected = set(['SNV', 'INDEL', 'CNA', 'SV'])
    if missing_text == '.':
        return expected
    missing = set(missing_text.split(','))
    return expected-missing

def get_handling(row: pd.Series) -> str:
    if row['primary'] == 0 and row['secondary'] == 0:
        return 'ignore'
    elif row['primary'] == 1 and row['secondary'] == 0:
        return 'single_primary'
    elif row['primary'] >= 2 and row['secondary'] == 0:
        return 'multi_primary'
    elif row['primary'] == 0 and row['secondary'] == 1:
        return 'single_met'
    elif row['primary'] == 0 and row['secondary'] >= 2:
        return 'multi_met'
    elif row['primary'] == 1 and row['secondary'] == 1:
        return 'single_primary_single_met'
    elif row['primary'] == 1 and row['secondary'] >= 2:
        return 'single_primary_multi_met'
    raise NotImplementedError(row)

def annotate_assignment_verdict_dynamic(table: pd.DataFrame, hframe: pd.DataFrame, verbose: bool=False) -> pd.DataFrame:
    dpc_weird_donors = set(['PPCG0388', 'PPCG0395', 'PPCG0425', 'PPCG0427', 'PPCG0428', 'PPCG0434', 'PPCG0442', 'PPCG1062', 'PPCG1072'])
    
    df = table.copy()
    handling_LUT = hframe.set_index('label')['handling'].to_dict()
    df['label'] = df['donor'] + '|' + df['vclass']
    df['handling'] = df['label'].map(handling_LUT)
    assert df['handling'].isna().sum() == 0

    mask0 = df['cohort']=='COMBI'
    mask1 = df['handling']=='single_primary_multi_met'
    mask2 = df['donor'].isin(dpc_weird_donors)
    
    mask_dpc = mask0 & mask1 & ~mask2
    mask_var = ~mask0 | ~mask1 | mask2 

    df.loc[mask_dpc, 'asmt'] = df.loc[mask_dpc, 'DPC_origin'] + ':' + df.loc[mask_dpc, 'DPC_ctype']
    df.loc[mask_var, 'asmt'] = df.loc[mask_var, 'var_origin'] + ':' + df.loc[mask_var, 'var_ctype']
    
    df.loc[mask_dpc, 'asmt_src'] = 'DPC'
    df.loc[mask_var, 'asmt_src'] = 'VAR'
    
    df = df.drop(['label', 'handling'], axis=1)
    assert df['asmt'].isna().sum() == 0
    assert df['asmt_src'].isna().sum() == 0
    return df.copy()

def annotate_assignment_verdict_var(table: pd.DataFrame) -> pd.DataFrame:
    df = table.copy()
    df['asmt'] = df['var_origin'] + ':' + df['var_ctype']
    df['asmt_src'] = 'VAR'
    assert df['asmt'].isna().sum() == 0
    assert df['asmt_src'].isna().sum() == 0
    return df.copy()

