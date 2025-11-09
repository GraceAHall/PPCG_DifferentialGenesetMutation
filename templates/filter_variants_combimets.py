
import argparse 
import pandas as pd 
import networkx as nx 
from typing import Tuple

from utils import filter_hypermutators

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
pd.set_option('display.max_colwidth', 40)
pd.options.display.float_format = '{:.2f}'.format


import warnings 
warnings.filterwarnings('ignore')


def main() -> None:
    args = load_cmdline_args()
    
    # load mutations 
    df = pd.read_csv(args.mutations, sep='\t', header=0)

    # drop duplicates 
    df = filter_duplicates(df)

    # stratify multi-sample donors
    counts = df.groupby('donor')['sample'].nunique()
    multi_sample_donors = counts[counts>=2].index.to_list()
    df_multi = df[df['donor'].isin(multi_sample_donors)]
    df_single = df[~df['donor'].isin(multi_sample_donors)]
    print(f"{len(multi_sample_donors)} donors have 2+ samples.")

    # assign variants to clones for multi sample donors
    cframe = assign_clones(df_multi, args.trees_dir, args.ccfs_dir, args.mettraj_clones, args.samplesheet)

    # label non- multi sample donor variants as all 'primary'
    df_single['label'] = 'primary'

    # label multi sample donors via transferring labels from sframe
    df_multi['ident'] = df_multi['donor'] + '|' + df_multi['vclass'] + '|' + df_multi['coords'] + '|' + df_multi['gene']
    cframe['ident']   = cframe['donor'] +   '|' + cframe['vclass'] +   '|' + cframe['coords'] +   '|' + cframe['gene']
    df_multi = df_multi.drop_duplicates(subset=['ident'])  # this removes sample-level duplicate records
    df_multi = df_multi.set_index('ident')
    cframe   = cframe.set_index('ident')
    df_multi['label'] = cframe['label'] # this does labelling.
    df_multi = df_multi.reset_index()
    df_multi = df_multi.drop('ident', axis=1)

    # merge 
    print('\nsingle sample labels:')
    print(df_single['label'].value_counts(dropna=False))
    print('\nmulti sample labels:')
    print(df_multi['label'].value_counts(dropna=False))
    print()
    print('WNT single: ', df_single[df_single['gene'].str.startswith('WNT')]['donor'].nunique())
    print('WNT multi: ', df_multi[df_multi['gene'].str.startswith('WNT')]['donor'].nunique())

    df = pd.concat([df_single, df_multi], ignore_index=True)
    df = df.drop('sample', axis=1)

    # actual filtering occurs here 
    df = df[df['label'].isin(['primary', 'mettraj'])].copy()
    df = df.drop('label', axis=1)
    df = df.sort_values(['donor', 'coords', 'gene'])

    # hypermutators
    df_filt, results = filter_hypermutators(df, args.zscore_thresh)
    results = results.reset_index()

    # write to file.
    df_filt.to_csv(args.outfile_mutations, sep='\t', index=False)
    results.to_csv(args.outfile_summary, sep='\t', index=False)


def filter_duplicates(table: pd.DataFrame) -> pd.DataFrame:
    # drop shallow deletions when the donor has a deep deletion
    df = table.copy()
    order_lut = {'CNA↓↓': 0, 'CNA↓': 1, 'CNA↑': 2}
    df['order'] = df['vtype'].apply(lambda x: order_lut.get(x, 10))
    df = df.sort_values('order')
    df = df.drop('order', axis=1)
    counts = pd.DataFrame(index=sorted(list(df['vclass'].unique())))
    counts['before'] = df['vclass'].value_counts()
    df = df.drop_duplicates(subset=['sample', 'coords', 'gene', 'vclass'])
    counts['after'] = df['vclass'].value_counts()
    counts['removed'] = counts['before'] - counts['after']
    print(counts)
    return df 



########################
### CLONE ASSIGNMENT ###
########################

class CloneAssigner:
    def __init__(self, T: nx.DiGraph, ccfs: pd.DataFrame, asmts: dict[str, str], purities: dict[str, float]) -> None:
        self.T = T
        self.ccfs = ccfs
        self.asmts = asmts
        self.purities = purities

    def assign(self, coords: str, obs: dict[str, float|None]) -> Tuple[str, str]:
        if coords in self.asmts:
            return self.assign_dpclust_direct(coords)
        return self.site_parsimony(obs)
    
    def assign_dpclust_direct(self, coords: str) -> Tuple[str, str]:
        clone = self.asmts[coords]
        meth = 'dpclust direct'
        return clone, meth

    def site_parsimony(self, obs: dict[str, float|None]) -> Tuple[str, str]:
        exact = True
        clones = self._parsimonious_clones_exact(set(obs.keys()))
        if len(clones) == 0:
            exact = False 
            clones = self._parsimonious_clones_inexact(set(obs.keys()))

        if exact and len(clones) >= 2:
            meth = 'site_parsimony_exact_multi'
        elif exact and len(clones) == 1:
            meth = 'site_parsimony_exact_single'
        elif len(clones) >= 2:
            meth = 'site_parsimony_inexact_multi'
        elif len(clones) == 1:
            meth = 'site_parsimony_inexact_single'
        else:
            print()
            print(obs)
            print(clones)
            raise RuntimeError

        if len(clones) == 1:
            return clones.pop(), meth
        return self._best_matching_clone(clones, obs), meth

    def _parsimonious_clones_exact(self, obs: set[str]) -> set[str]:
        clones = set()
        for node in nx.topological_sort(self.T):
            if obs == self.T.nodes[node]['samples']:
                clones.add(node)
        return clones

    def _parsimonious_clones_inexact(self, obs: set[str]) -> set[str]:

        def _trace(obs: set[str], leaf: str) -> str:
            node = leaf
            samples = self.T.nodes[node]['samples']
            while len(obs - samples) > 0:
                parents = list(self.T.predecessors(node))
                assert len(parents) == 1
                node = parents[0]
                samples = self.T.nodes[node]['samples']
            return node
        
        clones = set()
        leaves = [n for n in self.T.nodes() if self.T.out_degree(n)==0]
        for leaf in leaves:
            clones.add(_trace(obs, leaf))

        data = [(c, len(self.T.nodes[c]['samples'] - obs)) for c in clones]
        data = sorted(data, key=lambda x: x[1])
        min_extras = data[0][1]
        return set([x[0] for x in data if x[1]==min_extras])
    
    def _best_matching_clone(self, clones: set[str], obs: dict[str, float|None]) -> str:


        data = []
        for clone in clones:
            data.append((
                clone, 
                self._mean_ccf_diff(clone, obs),
                self._mean_unobs_purity(clone, obs)
            ))
        df = pd.DataFrame.from_records(data, columns=['clone', 'mean_ccf_diff', 'mean_unobs_purity'])
        df = df.set_index('clone')
        df['score'] = df.sum(axis=1)
        df = df.sort_values('score')
        
        # print()
        # print()
        # print(f"candidates: {clones}")
        # print(f"observed:   {obs}")
        # print(f"purity:     {self.purities}")
        # print()
        # print(self.ccfs)
        # print()
        # print(df)
        return df.index.to_list()[0]

    def _mean_unobs_purity(self, clone: str, obs: dict[str, float|None]) -> float:
        clone_samples = self.T.nodes[clone]['samples']
        obs_samples = set(obs.keys())
        unobs_samples = clone_samples - obs_samples
        if len(unobs_samples) == 0:
            return 0.0
        else:
            return sum([self.purities[s] for s in unobs_samples]) / len(unobs_samples)

    def _mean_ccf_diff(self, clone: str, obs: dict[str, float|None]) -> float:
        data = []
        for obs_sam, obs_ccf in obs.items():
            # comparing against SV clone
            if clone not in self.ccfs.index:
                data.append(0.5)

            # comparing SV against SNV clone
            elif obs_ccf is None:
                data.append(0.5)

            # comparing SNV/INDEL/CNA against SNV clone
            else:
                clust_ccf = self.ccfs.loc[clone, obs_sam]
                data.append(abs(clust_ccf-obs_ccf))
        return sum(data) / len(data)
    

def load_ccfs(filepath: str) -> pd.DataFrame:
    df = pd.read_csv(filepath, sep=',', header=0)
    df = df.rename(columns={'Cluster': 'cluster'})
    df['cluster'] = df['cluster'].apply(lambda x: str(x).replace('_', ''))
    df = df.set_index('cluster')
    df = df[[x for x in df.columns if x.startswith('PPCG')]]
    df.columns = [x.replace('_DNA', '') for x in df.columns]
    return df

def load_assignments(filepath: str) -> dict[str, str]:
    df = pd.read_csv(filepath, header=0)
    df.columns = [x.lower() if not x.startswith('PPCG') else x for x in df.columns]
    df.columns = [x.replace('_DNA', '') for x in df.columns]
    df['cluster'] = df['cluster'].astype(str)
    df['cluster'] = df['cluster'].apply(lambda x: x.replace('_', ''))
    df['chr'] = df['chr'].astype(str)
    df['pos'] = df['pos'].astype(str)
    df['coords'] = df['chr'] + ':' + df['pos']
    return df.set_index('coords')['cluster'].to_dict()

def load_phylo(filepath: str) -> pd.DataFrame:
    df = pd.read_csv(filepath, sep='\t', header=0)
    df['id'] = df['id'].apply(lambda x: str(x).replace('_', ''))
    df['parent'] = df['parent'].apply(lambda x: str(x).replace('_', ''))
    return df

def load_purities(filepath: str, donor: str) -> dict:
    df = pd.read_csv(filepath, sep='\t', header=0)
    df = df[df['donor']==donor].copy()
    return df.set_index('sample')['purity'].to_dict()

def load_clonetree(ccfs: pd.DataFrame, phylo: pd.DataFrame) -> nx.DiGraph:
    # init graph
    T = nx.DiGraph()

    # add all nodes and edges
    for rec in phylo.itertuples():
        T.add_edge(rec.parent, rec.id)

    # remove the root
    T.remove_node('0')

    # add observed samples for DPClust SNV clusters
    for node in T.nodes():
        if node not in ccfs.index:
            T.nodes[node]['samples'] = {}
            T.nodes[node]['ctype'] = 'SV'
            continue
        row = ccfs.loc[node]
        row = row[row>0]
        T.nodes[node]['samples'] = set(row.index.to_list())
        T.nodes[node]['ctype'] = 'SNV'
    
    # add observed samples SV-only clusters (reachable samples)
    for node in nx.topological_sort(T):
        if len(T.nodes[node]['samples']) == 0:
            reachable = set()
            for desc_node in nx.descendants(T, node):
                reachable.update(T.nodes[desc_node]['samples'])
            if len(reachable) == 0:
                parents = list(T.predecessors(node))
                assert len(parents) == 1
                assert len(T.nodes[parents[0]]['samples']) > 0
                reachable = T.nodes[parents[0]]['samples']
            T.nodes[node]['samples'] = reachable
    return T

def assign_clones(table: pd.DataFrame, trees_dir: str, ccfs_dir: str, mettraj_clones_path: str, samplesheet_path: str) -> pd.DataFrame:
    meta = pd.read_csv(mettraj_clones_path, sep='\t', header=0)
    meta = meta.set_index('donor')
    
    df = table.copy()
    df['est_ccf'] = df['est_ccf'].fillna(value='NA')

    data = []
    for donor in sorted(list(df['donor'].unique())):
        dfslice = df[df['donor']==donor].copy()
        met_clones = str(meta.loc[donor, 'met_trajectory_clones']).split(',')
        print(donor)

        # paths 
        phylo_path = f"{trees_dir}/{donor}.txt"
        ccfs_path = f"{ccfs_dir}/{donor}_Cluster_CCFs.csv"
        asmt_path = f"{ccfs_dir}/{donor}_SNV_CCF_Cluster_assignment.csv"
        outfile_tree = f"/home/grace/work/PPCGG_GeneSetPipeline/temp/angel_trees/{donor}.png"

        # supporting data
        ccfs = load_ccfs(ccfs_path)
        phylo = load_phylo(phylo_path)
        asmts = load_assignments(asmt_path)
        purities = load_purities(samplesheet_path, donor)

        # tree
        T = load_clonetree(ccfs, phylo)

        # assigner class
        assigner = CloneAssigner(T, ccfs, asmts, purities)

        for rec in dfslice.drop_duplicates(subset=['vclass', 'coords', 'gene']).itertuples():
            vclass = str(rec.vclass)
            coords = str(rec.coords)
            gene = str(rec.gene)
            varrecs = dfslice[
                (dfslice['vclass']==vclass) & 
                (dfslice['coords']==coords) & 
                (dfslice['gene']==gene)
            ].copy()

            obs = {}
            for rec in varrecs.itertuples():
                if rec.est_ccf == 'NA':
                    obs[rec.sample] = None
                else:
                    obs[rec.sample] = rec.est_ccf
            
            clone, method = assigner.assign(coords, obs)
            is_mettraj = 'mettraj' if clone in met_clones else 'other'
            data.append((donor, vclass, coords, gene, clone, method, is_mettraj))

    print('\n\n--- SUMMARY ---')
    cframe = pd.DataFrame.from_records(data, columns=['donor', 'vclass', 'coords', 'gene', 'clone', 'method', 'label'])
    temp = cframe.drop_duplicates(subset=['vclass', 'coords'])
    print()
    print(temp.groupby('vclass')['method'].value_counts().unstack().fillna(0).astype(int))
    print()
    print(temp.groupby('donor')['method'].value_counts().unstack().fillna(0).astype(int))
    print()
    print(temp.groupby('donor')['clone'].value_counts())
    print()
    print(temp.groupby('donor')['vclass'].value_counts().unstack().fillna(0).astype(int))
    print()
    pframe = temp.groupby('donor')['label'].value_counts(normalize=True).unstack().fillna(0.0).astype(float)
    # pframe
    pframe['phylo_mettraj'] = pd.NA 
    pframe['phylo_other'] = pd.NA
    for donor in sorted(list(df['donor'].unique())):
        met_clones = str(meta.loc[donor, 'met_trajectory_clones']).split(',')
        phylo_path = f"{trees_dir}/{donor}.txt"
        phylo = load_phylo(phylo_path)
        phylo = phylo.set_index('id')
        n_total = phylo['n_mut'].sum()
        n_mettraj = phylo.loc[met_clones, 'n_mut'].sum()
        prop_mettraj = n_mettraj / n_total
        prop_other = (n_total-n_mettraj) / n_total
        pframe.loc[donor, 'phylo_mettraj'] = prop_mettraj
        pframe.loc[donor, 'phylo_other'] = prop_other
    pframe = pframe.rename(columns={'mettraj': 'assigned_mettraj', 'other': 'assigned_other'})
    print(pframe)
    print()
    print(cframe.head())
    print()
    return cframe


def load_cmdline_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Merges all variant types into single file for a given PPCG donor.')
    parser.add_argument('--mutations', type=str, required=True, help='Path to merged variants (tsv).')
    parser.add_argument('--samplesheet', type=str, required=True, help='Path to samplesheet (used for purity information).')
    parser.add_argument('--trees-dir', type=str, required=True, help='Path to directory containing machina tree data.')
    parser.add_argument('--ccfs-dir', type=str, required=True, help='Path to directory containing dpclust clone ccfs.')
    parser.add_argument('--mettraj-clones', type=str, required=True, help='Path to file mapping donors to clones in metastatic trajectory.')
    parser.add_argument('--zscore-thresh', type=float, required=True, help='outlier threshold via Modified Z-Score (MAD) (standard deviations from the mean).')
    parser.add_argument('--outfile-mutations', type=str, required=True, help='Path to output .tsv file containing filtered variants.')
    parser.add_argument('--outfile-summary', type=str, required=True, help='Path to output .tsv file containing hypermutator summary.')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    main()

