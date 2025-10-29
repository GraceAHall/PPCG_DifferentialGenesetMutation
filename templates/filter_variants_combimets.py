
import argparse 
import pandas as pd 
import networkx as nx 
import seaborn as sns 

from typing import Tuple
from abc import ABC, abstractmethod
from collections import defaultdict
import matplotlib.pyplot as plt 

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

    # write to file. 
    df.to_csv(args.outfile, sep='\t', index=False)


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

        # met_clones = str(meta.loc[donor, 'met_trajectory_clones']).split(',')
        # oth_clones = [n for n in T.nodes() if n not in met_clones]
        # pos = nx.drawing.nx_pydot.graphviz_layout(T, prog='dot')
        # nx.draw_networkx_nodes(T, pos, nodelist=met_clones, node_color='salmon', node_size=1500)
        # nx.draw_networkx_nodes(T, pos, nodelist=oth_clones, node_color='skyblue', node_size=1500)
        # nx.draw_networkx_edges(T, pos, edge_color='gray')
        # nx.draw_networkx_labels(T, pos, labels={n: f"{n}\n{T.nodes[n]['ctype']}" for n in T.nodes()})
        # plt.savefig(outfile_tree)
        # plt.close()


def load_cmdline_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Merges all variant types into single file for a given PPCG donor.')
    parser.add_argument('--mutations', type=str, required=True, help='Path to merged variants (tsv).')
    parser.add_argument('--samplesheet', type=str, required=True, help='Path to samplesheet (used for purity information).')
    parser.add_argument('--trees-dir', type=str, required=True, help='Path to directory containing machina tree data.')
    parser.add_argument('--ccfs-dir', type=str, required=True, help='Path to directory containing dpclust clone ccfs.')
    parser.add_argument('--mettraj-clones', type=str, required=True, help='Path to file mapping donors to clones in metastatic trajectory.')
    parser.add_argument('--outfile', type=str, required=True, help='Path to output file containing filtered variants.')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    main()




# def load_clonetree(filepath: str) -> nx.DiGraph:
#     # need clone, samples


#     # df = pd.read_csv(f'/home/grace/work/PPCGG_GeneSetPipeline/data/phylogeny_angel/trees/{donor}.txt', sep='\t', header=0)
#     df = pd.read_csv(filepath, sep='\t', header=0)
#     df['id'] = df['id'].apply(lambda x: str(x).replace('_', ''))
#     df['parent'] = df['parent'].apply(lambda x: str(x).replace('_', ''))

#     T = nx.DiGraph()

#     # add nodes
#     T.add_node('0')
#     for rec in df.itertuples():
#         T.add_edge(rec.parent, rec.id)
#         if rec.id not in T.nodes():
#             T.add_node(rec.id, muts=rec.n_mut, status='unknown', type='sv', mettraj=False)

#     # add edges from phylo table
#     for rec in phylo.itertuples():
#         assert rec.id in T.nodes()
#         assert rec.parent in T.nodes()

#     T.remove_node('0')




# ##########################
# ### TREE INTROSPECTION ###
# ##########################

# # mark each variant as 'metastatic trajectory', 'private primary' or 'private secondary'. 
# # this uses the tree for the given donor. 
# def seeds(T: nx.DiGraph) -> set[Tuple[str, str]]:
#     # get all metastatic clones
    
#     def _is_valid_target(s1: str, s2: str) -> bool:
#         if s1 == s2:
#             return False
#         if s1 != 'prostate':
#             return False
#         return True
    
#     seedclones = set()
#     for n1, n2 in T.edges():
#         s1, s2 = T.nodes[n1]['site'], T.nodes[n2]['site']
#         c1, c2 = T.nodes[n1]['clone'], T.nodes[n2]['clone']
#         valid = _is_valid_target(s1, s2)
#         if valid and c1 == c2:
#             seedclones.add((c1, s2))
#     return seedclones

# def metastatic(T: nx.DiGraph, seeds: set[str]) -> set[str]:
#     if len(seeds) == 0:
#         return set()
    
#     # query primary clones to see whether they're along a path from root to filtered seedclones. 
#     pnodes = [n for n in nx.topological_sort(T) if T.nodes[n]['site']=='prostate']
#     root = pnodes[0]
    
#     trajectory = set()
#     for seedclone in seeds:
#         seednode = [n for n in pnodes if T.nodes[n]['clone']==seedclone][0]
#         path = nx.shortest_path(T, root, seednode)
#         for node in path:
#             clone = T.nodes[node]['clone']
#             trajectory.add(clone)
    
#     return trajectory

# def primary_private(T: nx.DiGraph, seeds: set[str]) -> set[str]:
#     metclones = metastatic(T, seeds)
#     pclones = set([T.nodes[n]['clone'] for n in nx.topological_sort(T) if T.nodes[n]['site']=='prostate'])
#     return pclones - metclones

# def secondary_private(T: nx.DiGraph, seeds: set[str]) -> dict[str, set[str]]:
#     info = defaultdict(set)
#     for n1, n2 in T.edges():
#         s1, s2 = T.nodes[n1]['site'], T.nodes[n2]['site']
#         c2 = T.nodes[n2]['clone']
#         if s1 == s2:
#             if s1 != 'prostate' and s2 != 'prostate':
#                 info[s2].add(c2)
#     return info


########################
### CLONE ASSIGNMENT ###
########################

# class CloneAssigner(ABC):

#     def __init__(self, table: pd.DataFrame, T: nx.DiGraph, ccfs: pd.DataFrame, asmts: pd.Series) -> None:
#         self.table = table
#         self.ccfs = ccfs
#         self.asmts = asmts
#         self.T = T
#         self.valid = set(self.ccfs.index.to_list())

#     @abstractmethod
#     def assign(self, coords: str, gene: str) -> str:
#         ...

#     def assign_site_parsimony_ccfs(self, df: pd.DataFrame) -> Tuple:
#         # get parsimonious clones (1+ candidates). clones will always correspond to a single sample. 
#         sample, candidate_clones, method = self.parsimonious_clones(set(df['sample'].unique()))
        
#         # get dpclust ccf of parimonious clones 
#         candidate_ccfs = [float(self.ccfs.loc[clone, sample]) for clone in candidate_clones]

#         # get estimated ccf of the mutation
#         obs_ccfs = df[df['sample']==sample]['est_ccf'].values
#         obs_ccf = 0 if len(obs_ccfs) != 1 else obs_ccfs[0]

#         # get the clone with most similar ccf to any candidate. 
#         clone, clone_ccf, _, diff = self._best_matching_clone_ccf(candidate_clones, candidate_ccfs, obs_ccf)
    
#         return clone, clone_ccf, obs_ccf, diff, method

#     def assign_site_parsimony_earliest(self, df: pd.DataFrame) -> Tuple:
#         # get parsimonious clones (1+ candidates). clones will always correspond to a single sample. 
#         sample, candidate_clones, method = self.parsimonious_clones(set(df['sample'].unique()))
#         for node in nx.topological_sort(self.T):
#             if node in candidate_clones:
#                 clone = self.T.nodes[node]['clone']
#                 obs_ccf = pd.NA
#                 clone_ccf = self.ccfs.loc[clone, sample]
#                 diff = abs(clone_ccf - obs_ccf)
#                 return clone, clone_ccf, obs_ccf, diff, method

#         raise RuntimeError
    
#     def parsimonious_clones(self, samples: set[str]) -> Tuple:
#         """
#         Assigns mutation to best fitting clone using observed pattern of anatomical sites. 
#         """
#         cloneset = self._parsimonious_clones_exact(samples)
#         if len(cloneset) >= 1:
#             c_clones = cloneset
#             method = 'site-parsimony exact'
#         else:
#             c_clones = self._parsimonious_clones_inexact(samples)
#             method = 'site-parsimony inexact'
        
#         # ensure 1+ candidate clones
#         c_samples = set([self.T.nodes[n]['sample'] for n in c_clones])
#         assert len(c_clones) >= 1
#         # ensure single sample for parsimonious clones 
#         if len(c_samples) != 1:
#             raise RuntimeError('2+ anatomical sites.')

#         sample = c_samples.pop()
#         candidate_clones = sorted(list(c_clones))
#         return sample, candidate_clones, method 

#     def _parsimonious_clones_exact(self, observed: set[str]) -> set[str]:
#         # correct method
#         T = self.T
#         cloneset = set()
#         for node in nx.topological_sort(T):
#             clone = T.nodes[node]['clone']
#             if clone not in self.valid:
#                 continue
            
#             # which samples are reachable from this clone? 
#             reachable = set()
#             reachable = reachable | set([T.nodes[n]['sample'] for n in T.nodes() if T.nodes[n]['clone']==clone]) # this clone's sample
#             reachable = reachable | set(T.nodes[n]['sample'] for n in nx.descendants(T, node)) # descendent clones samples
            
#             # if there is parsimony, add to candidates.
#             # ie the pattern of reachable samples matches the observed pattern of samples with the mutation. 
#             if reachable == observed:
#                 cloneset.add(clone)
#         return cloneset

#     def _parsimonious_clones_inexact(self, observed: set[str]) -> set[str]:
#         # shittier method for donors where MACHINA migration graph may have some inconsistency with the data. 
#         T = self.T
#         cloneset = set()

#         # query each node, recording reachable sites. 
#         # missing: sites which variant WAS observed in, but UNREACHABLE by the clone. 
#         # extra: sites which variant WAS NOT observed in, but REACHABLE by the clone. 
#         candidates = []
#         for node in sorted(list(nx.topological_sort(T))):
#             clone = T.nodes[node]['clone']
#             if clone not in self.valid:
#                 continue
#             if clone in candidates:
#                 continue 
#             reachable = set()
#             reachable = reachable | set([T.nodes[n]['sample'] for n in T.nodes() if T.nodes[n]['clone']==clone]) # this clone's sample
#             reachable = reachable | set(T.nodes[n]['sample'] for n in nx.descendants(T, node)) # descendent clones samples
#             missing = observed - reachable
#             extras = reachable - observed
            
#             # Assume pidgeonhole. 
#             # All observed samples must be reachable by clone. 
#             # Extras may be simply due to low variant CCF in some samples.  
#             # Especially true when the clone is in the primary prostate sample (and not truncal).
#             if len(missing) == 0:
#                 candidates.append((clone, extras))

#         # select candidate clones as those with minimal number of extra reachable samples. 
#         min_extras = min([len(x[1]) for x in candidates])  
#         for clone, extrasamples in candidates:
#             if len(extrasamples) == min_extras:
#                 cloneset.add(clone)
        
#         # Final choice: lowest candidate clone in the tree. 
#         # This is because CCF gets lower as we descend the tree. 
#         # Assumption is (as above), the issue was caused by low clone CCF in the primary prostate. 
#         # This is why we don't 'observe' the variant in the prostate. 
#         for node in sorted(list(nx.topological_sort(T)), reverse=True):
#             clone = T.nodes[node]['clone']
#             if clone in cloneset:
#                 return {clone}

#         raise RuntimeError
    
#     def _best_matching_clone_ccf(self, candidate_clones: list[str], candidate_ccfs: list[float], est_ccf: float) -> list:
#         data = []
#         for clone, c_ccf in zip(candidate_clones, candidate_ccfs):
#             data.append((clone, c_ccf, est_ccf))
#         df = pd.DataFrame.from_records(data=data, columns=['clone', 'clone_ccf', 'est_ccf'])
#         df['diff'] = df['clone_ccf'] - df['est_ccf']
#         df['diff'] = df['diff'].apply(lambda x: abs(x))
#         df = df.sort_values('diff')
#         return df.iloc[0].to_list()

# class SNVassigner(CloneAssigner):

#     def assign(self, coords: str, gene: str) -> Tuple:
#         df = self.table[(self.table['coords']==coords) & (self.table['gene']==gene)].copy()
#         info = self.assign_dpclust_direct(df, coords)
#         if info is not None:
#             return info
#         return self.assign_site_parsimony_ccfs(df)
    
#     def assign_dpclust_direct(self, df: pd.DataFrame, coords: str) -> Tuple|None:
#         if coords not in self.asmts:
#             return None
        
#         clone = self.asmts[coords]
#         sample = [self.T.nodes[n]['sample'] for n in self.T.nodes() if self.T.nodes[n]['clone']==clone][0]
#         obs_ccfs = df[df['sample']==sample]['est_ccf'].values
#         if len(obs_ccfs) == 0:
#             return None 
#         if len(obs_ccfs) >= 2:
#             raise RuntimeError

#         clone_ccf = self.ccfs.loc[clone, sample]
#         obs_ccf = obs_ccfs[0]
#         diff = abs(clone_ccf - obs_ccf)
#         method = 'dpclust direct'
#         return clone, clone_ccf, obs_ccf, diff, method
    
# class INDELassigner(CloneAssigner):

#     def assign(self, coords: str, gene: str) -> Tuple:
#         df = self.table[(self.table['coords']==coords) & (self.table['gene']==gene)].copy()
#         return self.assign_site_parsimony_ccfs(df)

# class SVassigner(CloneAssigner):

#     def assign(self, coords: str, gene: str) -> Tuple:
#         df = self.table[(self.table['coords']==coords) & (self.table['gene']==gene)].copy()
#         return self.assign_site_parsimony_earliest(df)

# class CNAassigner(CloneAssigner):

#     def assign(self, coords: str, gene: str) -> Tuple:
#         df = self.table[(self.table['coords']==coords) & (self.table['gene']==gene)].copy()
#         return self.assign_site_parsimony_ccfs(df)

