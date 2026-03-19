
import pandas as pd
from glob import glob 

##############
### DONORS ###
##############

def filter_ppcg_control(table: pd.DataFrame, verbose: bool=False) -> pd.DataFrame:
    valid = ['PPCG0001', 'PPCG0003', 'PPCG0004', 'PPCG0006', 'PPCG0008', 'PPCG0010', 'PPCG0017', 'PPCG0019', 'PPCG0022', 'PPCG0052', 'PPCG0053', 'PPCG0074', 'PPCG0134', 'PPCG0138', 'PPCG0159', 'PPCG0161', 'PPCG0172', 'PPCG0174', 'PPCG0203', 'PPCG0261', 'PPCG0264', 'PPCG0275', 'PPCG0280', 'PPCG0283', 'PPCG0284', 'PPCG0289', 'PPCG0306', 'PPCG0307', 'PPCG0314', 'PPCG0342', 'PPCG0347', 'PPCG0348', 'PPCG0349', 'PPCG0408', 'PPCG0443', 'PPCG0452', 'PPCG0487', 'PPCG0516', 'PPCG0517', 'PPCG0519', 'PPCG0525', 'PPCG0548', 'PPCG0553', 'PPCG0567', 'PPCG0568', 'PPCG0592', 'PPCG0594', 'PPCG0595', 'PPCG0612', 'PPCG0615', 'PPCG0631', 'PPCG0640', 'PPCG0696', 'PPCG0775', 'PPCG0792', 'PPCG0805', 'PPCG0813', 'PPCG0830', 'PPCG0831', 'PPCG0837', 'PPCG0885', 'PPCG0889', 'PPCG0910', 'PPCG0963', 'PPCG0977']
    df = table.copy()
    before = set(df['donor'].unique())
    df = df[df['donor'].isin(valid)]
    after = set(df['donor'].unique())
    if verbose:
        print('\nRemove donors not in the 65-donor control group')
        print(len(after), f"(removed {', '.join(sorted(list(before-after)))})")
    return df.copy()

def filter_donors_without_prostate(table: pd.DataFrame, verbose: bool=False) -> pd.DataFrame:
    # filter donors without any mutations observed in their primary tissue sample
    df = table.copy()
    before = set(df['donor'].unique())
    valid = set(df[df['tissue']=='Primary']['donor'].unique())
    df = df[df['donor'].isin(valid)]
    after = set(df['donor'].unique())
    if verbose:
        print('\nRemove donors without primary tissue')
        print(len(after), f"(removed {', '.join(sorted(list(before-after)))})")
    return df.copy()

def filter_donors_without_dpclust(table: pd.DataFrame, dpclust_dirpath: str, verbose: bool=False) -> pd.DataFrame:
    df = table.copy()
    
    # filter donors without DPClust CCF file. 
    before = set(df['donor'].unique())
    donor2ccfpath = {x.split('/')[-1][:8]: x for x in glob(f"{dpclust_dirpath}/*_Cluster_CCFs.csv")}
    valid = set(donor2ccfpath.keys())
    df = df[df['donor'].isin(valid)]
    after = set(df['donor'].unique())
    if verbose:
        print('\nRemove donors without DPClust')
        print(len(after), f"(removed {', '.join(sorted(list(before-after)))})")

    # filter donors missing a sample in their DPClust CCF file.
    before = after
    valid = set()
    for donor in df['donor'].unique():
        fields = pd.read_csv(donor2ccfpath[donor], sep=',', header=0).columns.to_list()
        dpc_samples = set([f[:9] for f in fields if f.startswith('PPCG')])
        mut_samples = set(df[df['donor']==donor]['sample'].unique())
        if len(mut_samples - dpc_samples) == 0:
            valid.add(donor)
    df = df[df['donor'].isin(valid)]
    after = set(df['donor'].unique())
    if verbose:
        print('\nRemove donors where DPClust is missing a sample')
        print(len(after), f"(removed {', '.join(sorted(list(before-after)))})")
    
    # filter donors which don't have a DPClust 'Cluster_Type' marked 'Trunk'.
    before = after
    valid = set()
    for donor in df['donor'].unique():
        ccfs = pd.read_csv(donor2ccfpath[donor], sep=',', header=0)
        ctypes = set(ccfs['Cluster_Type'].unique())
        if 'Trunk' in ctypes:
            valid.add(donor)
    df = df[df['donor'].isin(valid)]
    after = set(df['donor'].unique())
    if verbose:
        print('\nRemove donors where DPClust has no truncal cluster')
        print(len(after), f"(removed {', '.join(sorted(list(before-after)))})")

    return df.copy()

def filter_donors_without_trees(table: pd.DataFrame, conipher_dirpath: str, verbose: bool=False) -> pd.DataFrame:
    df = table.copy()
    before = set(df['donor'].unique())
    valid = set([x.split('/')[-2][:8] for x in glob(f"{conipher_dirpath}/*/allTrees.txt")])
    df = df[df['donor'].isin(valid)]
    after = set(df['donor'].unique())
    if verbose:
        print('\nRemove donors without CONIPHER tree')
        print(len(after), f"(removed {', '.join(sorted(list(before-after)))})")
    return df.copy()


#####################
### HYPERMUTATORS ###
#####################

def filter_hypermutators(table: pd.DataFrame, max_mutated_genes: int=500, verbose: bool=False) -> pd.DataFrame:
    df = table.copy()

    blacklist = set()
    # for vclass in ['SNV','INDEL','CNA','SV']:
    #     counts = df[df['vclass']==vclass].groupby('donor')['gene'].nunique()
    #     invalid = set(counts[counts>max_mutated_genes].index.to_list())
    #     blacklist.update(invalid)
    counts = df.groupby('donor')['gene'].nunique()
    invalid = set(counts[counts>max_mutated_genes].index.to_list())
    blacklist.update(invalid)
    
    before = set(df['donor'].unique())
    df = df[~df['donor'].isin(blacklist)]
    after = set(df['donor'].unique())
    
    if verbose:
        print(f'\nRemove hypermutators (more than {max_mutated_genes} mutated genes)')
        print(len(after), f"(removed {', '.join(sorted(list(before-after)))})")
    
    return df.copy()


#################
### MUTATIONS ###
#################

def filter_primary_leaves(table: pd.DataFrame, verbose: bool=False) -> pd.DataFrame:
    blacklist = ['Primary:Leaf']

    df = table.copy()
    before = set(df['ID'].unique())
    df = df[~df['asmt'].isin(blacklist)]
    after = set(df['ID'].unique())
    if verbose:
        print('\nRemove primary leaf clones')
        print(len(after), f"(removed {len(before-after)} mutations)")

    return df.copy()

def filter_secondary_leaves(table: pd.DataFrame, verbose: bool=False) -> pd.DataFrame:
    blacklist = ['Metastasis:Leaf']

    df = table.copy()
    before = set(df['ID'].unique())
    df = df[~df['asmt'].isin(blacklist)]
    after = set(df['ID'].unique())
    if verbose:
        print('\nRemove primary/secondary leaf clones')
        print(len(after), f"(removed {len(before-after)} mutations)")

    return df.copy()


