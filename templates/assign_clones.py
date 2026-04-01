
import argparse 
import warnings 
import pandas as pd
from glob import glob 

import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from site_parsimony import SiteParsimonyAssigner


warnings.filterwarnings('ignore')

pd.options.display.float_format = "{:,.3f}".format
pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', 200)
pd.set_option('display.width', 200)

# INFILE_MUTATIONS = '/home/grace/work/PPCG_DifferentialGenesetMutation/outputs/alldonors_01042026/variant_processing/mutations.filtered.tsv'
# INFILE_CLONE_META = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/angel/26.03.2026/metastatic_clones_with_ancestors.tsv'
# INDIR_DPCLUST = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/angel/26.03.2026/ccf/Clonal_trees'
# INDIR_CONIPHER = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/angel/06.03.2026/conipher_trees'
# OUTFILE_MUTATIONS = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/manual/mutations.assigned.010426.tsv'

def main():
    args = load_cmdline_args()
    muts = load_mutations(args.mutations)

    # split donors into non-assignment and assignment shards
    meta = pd.read_csv(args.clone_meta, sep='\t', header=0)
    assign_donors = set(meta['patient'].unique())

    # handle donors to skip
    shard1 = muts[~muts['donor'].isin(assign_donors)].copy()  # ignored 
    shard1['clone'] = pd.NA 

    # handle donors to assign
    shard2 = muts[muts['donor'].isin(assign_donors)].copy()   # to be assigned
    shard2 = assign_clones(shard2, args)

    # recombine
    asmts = pd.concat([shard1, shard2], ignore_index=True)
    asmts = asmts.sort_values(['donor', 'sample', 'gene'])
    asmts.drop('donor', axis=1).to_csv(args.outfile, sep='\t', index=False, float_format='%.2f')
    asmts.head(10)

def load_mutations(filepath: str) -> pd.DataFrame:
    df = pd.read_csv(filepath, sep='\t', header=0)
    return df 


def assign_clones(muts: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    # some samples are banned if don't appear in DPC cluster CCFs
    BANNED = ['PPCG0388a']
    df = muts[~muts['sample'].isin(BANNED)].copy()
    n_donors = df['donor'].nunique()

    i = 0
    table = pd.DataFrame()
    for donor, donor_df in df.groupby('donor'):
        print(f"processed {i}/{n_donors} donors...", end='\r')
        i += 1
        # print()
        # print(donor)
        tree_path = f"{args.trees_dir}/{donor}_conipher_tree/allTrees.txt"
        ccfs_path = f"{args.ccfs_dir}/{donor}_Cluster_CCFs.csv"
        asmt_path = f"{args.ccfs_dir}/{donor}_SNV_CCF_Cluster_assignment.csv"

        snvs = donor_df[donor_df['vclass']=='SNV'].copy()
        indels = donor_df[donor_df['vclass']=='INDEL'].copy()
        cna = donor_df[donor_df['vclass']=='CNA'].copy()
        svs = donor_df[donor_df['vclass']=='SV'].copy()

        snvs = assign_snvs(snvs, ccfs_path, tree_path, asmt_path)
        indels = assign_indels(indels, ccfs_path, tree_path)
        cna = assign_cna(cna, ccfs_path, tree_path)
        svs = assign_svs(svs, ccfs_path, tree_path)
        
        donor_df = pd.concat([snvs, indels, cna, svs], ignore_index=True)
        assert donor_df['clone'].isna().sum() == 0
        assert donor_df['method'].isna().sum() == 0
        # print(donor_df.groupby('vclass')['clone'].value_counts().unstack().fillna(0).astype(int))
        table = pd.concat([table, donor_df], ignore_index=True)
    print(f"processed {i}/{n_donors} donors... done.")
    assert table['clone'].isna().sum() == 0
    assert table['method'].isna().sum() == 0
    return table 


def assign_snvs(table: pd.DataFrame, ccfs_path: str, tree_path: str, asmt_path: str) -> pd.DataFrame:
    assert set(table['vclass'].unique()) == set(['SNV'])
    assert table['est_ccf'].isna().sum() == 0 

    def load_dpclust_asmts(filepath: str) -> dict:
        df = pd.read_csv(filepath, sep=',', header=0)
        df = df.dropna(subset=['Cluster'])
        df['Cluster'] = df['Cluster'].apply(lambda x: str(x).replace('_', '').split('.')[0])
        
        # hacky hot-fix for PPCG0435 where clone 11 was merged into clone 4
        donor = filepath.split('/')[-1][:8]
        if donor == 'PPCG0435':
            df['Cluster'] = df['Cluster'].replace('11', '4')

        df['chr'] = df['chr'].astype(str)
        df['pos'] = df['pos'].astype(str)
        df['ident'] = df['chr'] + ':' + df['pos']
        return df.set_index('ident')['Cluster'].to_dict()

    df = table.copy()

    # map DPC assigned snv chrpos to clusters, prepare mutations 
    chrpos2clust = load_dpclust_asmts(asmt_path)
    df['chrpos'] = df['coords'].apply(lambda x: x.split('-')[0])

    # no funny business
    if any(['chr' in x for x in chrpos2clust.keys()]):
        raise NotImplementedError
    if any(['chr' in x for x in df['chrpos'].unique()]):
        raise NotImplementedError

    # do direct DPC clone assignment 
    df['clone'] = df['chrpos'].map(chrpos2clust)
    
    # handle DPC direct assigned shard
    shard1 = df[df['clone'].notna()].copy()
    shard1['method'] = 'DPC direct'

    # handle DPC direct non-assigned shard
    shard2 = df[df['clone'].isna()].copy()
    assigner = SiteParsimonyAssigner(ccfs_path=ccfs_path, tree_path=tree_path)
    coords2clone = dict()
    coords2label = dict()
    for coords, coords_df in shard2.groupby('coords'):
        # because PPCG0332|6:139225153 (same location occupied by two genes). 
        # don't know if ccf would actually ever be different but whatever. 
        if coords_df.groupby('sample')['est_ccf'].nunique().max() >= 2:
            print(coords_df)
            raise NotImplementedError
        obs_df = coords_df.sort_values('est_ccf', ascending=False).drop_duplicates('sample')
        obs_dict = obs_df.set_index('sample')['est_ccf'].to_dict()
        clone, score, is_exact = assigner.assign(obs_samples2ccfs=obs_dict) # type: ignore
        coords2clone[coords] = clone
        coords2label[coords] = 'site-parsimony (exact)' if is_exact else 'site-parsimony (inexact)'
    shard2['clone'] = shard2['coords'].map(coords2clone)
    shard2['method'] = shard2['coords'].map(coords2label)

    # rejoin and return 
    df = pd.concat([shard1, shard2], ignore_index=False)
    df = df.drop('chrpos', axis=1)
    return df 

def assign_indels(table: pd.DataFrame, ccfs_path: str, tree_path: str) -> pd.DataFrame:
    assert set(table['vclass'].unique()) == set(['INDEL'])
    assert table['est_ccf'].isna().sum() == 0 
    
    df = table.copy()
    assigner = SiteParsimonyAssigner(ccfs_path=ccfs_path, tree_path=tree_path)
    coords2clone = dict()
    coords2label = dict()
    for coords, coords_df in df.groupby('coords'):
        obs_df = coords_df.sort_values('est_ccf', ascending=False).drop_duplicates('sample')
        obs_dict = obs_df.set_index('sample')['est_ccf'].to_dict()
        clone, score, is_exact = assigner.assign(obs_samples2ccfs=obs_dict) # type: ignore
        coords2clone[coords] = clone
        coords2label[coords] = 'site-parsimony (exact)' if is_exact else 'site-parsimony (inexact)'
    df['clone'] = df['coords'].map(coords2clone)
    df['method'] = df['coords'].map(coords2label)
    return df

def assign_svs(table: pd.DataFrame, ccfs_path: str, tree_path: str) -> pd.DataFrame:
    assert set(table['vclass'].unique()) == set(['SV'])
    assert table['est_ccf'].notna().sum() == 0 

    df = table.copy()
    assigner = SiteParsimonyAssigner(ccfs_path=ccfs_path, tree_path=tree_path)
    gene2clone = dict()
    gene2label = dict()
    for gene, gene_df in df.groupby('gene'):
        obs_df = gene_df.drop_duplicates('sample').copy()
        obs_df['est_ccf'] = 1.0
        obs_dict = obs_df.set_index('sample')['est_ccf'].to_dict()
        clone, score, is_exact = assigner.assign(obs_samples2ccfs=obs_dict) # type: ignore
        gene2clone[gene] = clone
        gene2label[gene] = 'site-parsimony (exact)' if is_exact else 'site-parsimony (inexact)'
    df['clone'] = df['gene'].map(gene2clone)
    df['method'] = df['gene'].map(gene2label)

    return df

def assign_cna(table: pd.DataFrame, ccfs_path: str, tree_path: str) -> pd.DataFrame:
    assert set(table['vclass'].unique()) == set(['CNA'])
    assert table['est_ccf'].isna().sum() == 0 

    df = table.copy()
    assigner = SiteParsimonyAssigner(ccfs_path=ccfs_path, tree_path=tree_path)
    gene2clone = dict()
    gene2label = dict()
    for gene, gene_df in df.groupby('gene'):
        obs_df = gene_df.sort_values('est_ccf', ascending=False).drop_duplicates('sample')
        obs_dict = obs_df.set_index('sample')['est_ccf'].to_dict()
        clone, score, is_exact = assigner.assign(obs_samples2ccfs=obs_dict) # type: ignore
        gene2clone[gene] = clone
        gene2label[gene] = 'site-parsimony (exact)' if is_exact else 'site-parsimony (inexact)'
    df['clone'] = df['gene'].map(gene2clone)
    df['method'] = df['gene'].map(gene2label)
    return df



def load_cmdline_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Merges all variant types into single file for a given PPCG donor.')
    parser.add_argument('--mutations', type=str, required=True, help='Path to merged variants (tsv).')
    parser.add_argument('--clone-meta', type=str, required=True, help='Path to samplesheet (used for purity information).')
    parser.add_argument('--trees-dir', type=str, required=True, help='Path to directory containing CONIPHER trees.')
    parser.add_argument('--ccfs-dir', type=str, required=True, help='Path to directory containing DPClust ccfs.')
    parser.add_argument('--outfile', type=str, required=True, help='Path to output .tsv file containing filtered variants.')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    main()

