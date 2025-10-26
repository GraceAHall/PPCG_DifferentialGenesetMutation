
import argparse 
import pandas as pd 
import networkx as nx 
import seaborn as sns 

from typing import Tuple
from abc import ABC, abstractmethod
from collections import defaultdict

from utils import load_clonetree, modify_structure, load_dpclust_ccfs, load_dpclust_assignments

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
    cframe = assign_clones(df_multi, args.trees_dir, args.clones_dir)

    # assign labels for each clone assignment
    cframe = assign_labels(cframe, args.trees_dir)

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
    # print('WNT single: ', df_single[df_single['gene'].str.startswith('WNT')]['donor'].nunique())
    # print('WNT multi: ', df_multi[df_multi['gene'].str.startswith('WNT')]['donor'].nunique())

    df = pd.concat([df_single, df_multi], ignore_index=True)
    df = df.drop('sample', axis=1)

    # actual filtering occurs here 
    df = df[df['label'].isin(['primary', 'metastatic trajectory'])].copy()
    df = df.drop('label', axis=1)

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

def assign_clones(df: pd.DataFrame, trees_dir: str, clones_dir: str) -> pd.DataFrame:
    assigners_l = [
        ('SNV', SNVassigner),
        ('INDEL', INDELassigner),
        ('SV', SVassigner),
        ('CNA', CNAassigner),
    ]
    data = []
    for donor in sorted(list(df['donor'].unique())):
        # paths 
        tree_labels = f"{trees_dir}/clone_labels/{donor}.labeling"
        tree_edges = f"{trees_dir}/clone_edges/{donor}.tree"
        tree_dotfile = f"{trees_dir}/clone_dots/{donor}.dot"
        sam2site = f"{trees_dir}/sample2site/{donor}.tsv"
        clone_ccfs = f"{clones_dir}/{donor}_Cluster_CCFs.csv"
        clone_asmt = f"{clones_dir}/{donor}_SNV_CCF_Cluster_assignment.csv"
        
        # load clonetree and DPclust clone ccfs 
        T = load_clonetree(tree_labels, tree_edges, tree_dotfile, sam2site)
        T = modify_structure(T)
        ccfs = load_dpclust_ccfs(clone_ccfs)
        asmts = load_dpclust_assignments(clone_asmt)
        
        # assign variants of each class to clones
        for vclass, assigner_c in assigners_l:
            dfslice = df[(df['donor']==donor) & (df['vclass']==vclass)].copy()
            assigner = assigner_c(dfslice, T, ccfs, asmts)
            for rec in dfslice.drop_duplicates(subset=['coords', 'gene']).itertuples():
                coords = str(rec.coords)
                gene = str(rec.gene)
                clone, clone_ccf, est_ccf, diff, method = assigner.assign(coords, gene)
                data.append((donor, vclass, coords, gene, clone, clone_ccf, est_ccf, diff, method))

    cframe = pd.DataFrame.from_records(data, columns=['donor', 'vclass', 'coords', 'gene', 'clone', 'clone_ccf', 'est_ccf', 'diff', 'method'])
    print(cframe.groupby('vclass')['method'].value_counts())
    print()
    print(cframe.head())
    return cframe

def assign_labels(df: pd.DataFrame, trees_dir: str) -> pd.DataFrame:
    lut = defaultdict(dict)
    for donor in sorted(list(df['donor'].unique())):
        # load clonetree 
        tree_labels = f"{trees_dir}/clone_labels/{donor}.labeling"
        tree_edges = f"{trees_dir}/clone_edges/{donor}.tree"
        tree_dotfile = f"{trees_dir}/clone_dots/{donor}.dot"
        sam2site = f"{trees_dir}/sample2site/{donor}.tsv"
        T = load_clonetree(tree_labels, tree_edges, tree_dotfile, sam2site)
        T = modify_structure(T)

        # create lut mapping clones to labels 
        seed_info = seeds(T)
        seed_clones = set([x[0] for x in seed_info])
        mt_clones = metastatic(T, seed_clones)
        pp_clones = primary_private(T, seed_clones)
        sp_info = secondary_private(T, seed_clones)
        sp_clones = set()
        for site, clones in sp_info.items():
            sp_clones.update(list(clones))
        assert len(mt_clones & pp_clones) == 0
        assert len(mt_clones & sp_clones) == 0
        assert len(pp_clones & sp_clones) == 0
        for clone in mt_clones:
            lut[donor][clone] = 'metastatic trajectory'
        for clone in pp_clones:
            lut[donor][clone] = 'private primary'
        for clone in sp_clones:
            lut[donor][clone] = 'private secondary'
        
        assert len(lut[donor]) == len(set([T.nodes[n]['clone'] for n in T.nodes()]))

    df['clone'] = df['clone'].astype(str)
    for idx, row in df.iterrows():
        df.loc[idx, 'label'] = lut[row['donor']][row['clone']]

    print(df['clone'].isna().sum())
    print(df['label'].isna().sum())
    print(df.head())
    return df 



##########################
### TREE INTROSPECTION ###
##########################

# mark each variant as 'metastatic trajectory', 'private primary' or 'private secondary'. 
# this uses the tree for the given donor. 
def seeds(T: nx.DiGraph) -> set[Tuple[str, str]]:
    # get all metastatic clones
    
    def _is_valid_target(s1: str, s2: str) -> bool:
        if s1 == s2:
            return False
        if s1 != 'prostate':
            return False
        return True
    
    seedclones = set()
    for n1, n2 in T.edges():
        s1, s2 = T.nodes[n1]['site'], T.nodes[n2]['site']
        c1, c2 = T.nodes[n1]['clone'], T.nodes[n2]['clone']
        valid = _is_valid_target(s1, s2)
        if valid and c1 == c2:
            seedclones.add((c1, s2))
    return seedclones

def metastatic(T: nx.DiGraph, seeds: set[str]) -> set[str]:
    if len(seeds) == 0:
        return set()
    
    # query primary clones to see whether they're along a path from root to filtered seedclones. 
    pnodes = [n for n in nx.topological_sort(T) if T.nodes[n]['site']=='prostate']
    root = pnodes[0]
    
    trajectory = set()
    for seedclone in seeds:
        seednode = [n for n in pnodes if T.nodes[n]['clone']==seedclone][0]
        path = nx.shortest_path(T, root, seednode)
        for node in path:
            clone = T.nodes[node]['clone']
            trajectory.add(clone)
    
    return trajectory

def primary_private(T: nx.DiGraph, seeds: set[str]) -> set[str]:
    metclones = metastatic(T, seeds)
    pclones = set([T.nodes[n]['clone'] for n in nx.topological_sort(T) if T.nodes[n]['site']=='prostate'])
    return pclones - metclones

def secondary_private(T: nx.DiGraph, seeds: set[str]) -> dict[str, set[str]]:
    info = defaultdict(set)
    for n1, n2 in T.edges():
        s1, s2 = T.nodes[n1]['site'], T.nodes[n2]['site']
        c2 = T.nodes[n2]['clone']
        if s1 == s2:
            if s1 != 'prostate' and s2 != 'prostate':
                info[s2].add(c2)
    return info


########################
### CLONE ASSIGNMENT ###
########################

class CloneAssigner(ABC):

    def __init__(self, table: pd.DataFrame, T: nx.DiGraph, ccfs: pd.DataFrame, asmts: pd.Series) -> None:
        self.table = table
        self.ccfs = ccfs
        self.asmts = asmts
        self.T = T
        self.valid = set(self.ccfs.index.to_list())

    @abstractmethod
    def assign(self, coords: str, gene: str) -> str:
        ...

    def assign_site_parsimony_ccfs(self, df: pd.DataFrame) -> Tuple:
        # get parsimonious clones (1+ candidates). clones will always correspond to a single sample. 
        sample, candidate_clones, method = self.parsimonious_clones(set(df['sample'].unique()))
        
        # get dpclust ccf of parimonious clones 
        candidate_ccfs = [float(self.ccfs.loc[clone, sample]) for clone in candidate_clones]

        # get estimated ccf of the mutation
        obs_ccfs = df[df['sample']==sample]['est_ccf'].values
        obs_ccf = 0 if len(obs_ccfs) != 1 else obs_ccfs[0]

        # get the clone with most similar ccf to any candidate. 
        clone, clone_ccf, _, diff = self._best_matching_clone_ccf(candidate_clones, candidate_ccfs, obs_ccf)
    
        return clone, clone_ccf, obs_ccf, diff, method

    def assign_site_parsimony_earliest(self, df: pd.DataFrame) -> Tuple:
        # get parsimonious clones (1+ candidates). clones will always correspond to a single sample. 
        sample, candidate_clones, method = self.parsimonious_clones(set(df['sample'].unique()))
        for node in nx.topological_sort(self.T):
            if node in candidate_clones:
                clone = self.T.nodes[node]['clone']
                obs_ccf = pd.NA
                clone_ccf = self.ccfs.loc[clone, sample]
                diff = abs(clone_ccf - obs_ccf)
                return clone, clone_ccf, obs_ccf, diff, method

        raise RuntimeError
    
    def parsimonious_clones(self, samples: set[str]) -> Tuple:
        """
        Assigns mutation to best fitting clone using observed pattern of anatomical sites. 
        """
        cloneset = self._parsimonious_clones_exact(samples)
        if len(cloneset) >= 1:
            c_clones = cloneset
            method = 'site-parsimony exact'
        else:
            c_clones = self._parsimonious_clones_inexact(samples)
            method = 'site-parsimony inexact'
        
        # ensure 1+ candidate clones
        c_samples = set([self.T.nodes[n]['sample'] for n in c_clones])
        assert len(c_clones) >= 1
        # ensure single sample for parsimonious clones 
        if len(c_samples) != 1:
            raise RuntimeError('2+ anatomical sites.')

        sample = c_samples.pop()
        candidate_clones = sorted(list(c_clones))
        return sample, candidate_clones, method 

    def _parsimonious_clones_exact(self, observed: set[str]) -> set[str]:
        # correct method
        T = self.T
        cloneset = set()
        for node in nx.topological_sort(T):
            clone = T.nodes[node]['clone']
            if clone not in self.valid:
                continue
            
            # which samples are reachable from this clone? 
            reachable = set()
            reachable = reachable | set([T.nodes[n]['sample'] for n in T.nodes() if T.nodes[n]['clone']==clone]) # this clone's sample
            reachable = reachable | set(T.nodes[n]['sample'] for n in nx.descendants(T, node)) # descendent clones samples
            
            # if there is parsimony, add to candidates.
            # ie the pattern of reachable samples matches the observed pattern of samples with the mutation. 
            if reachable == observed:
                cloneset.add(clone)
        return cloneset

    def _parsimonious_clones_inexact(self, observed: set[str]) -> set[str]:
        # shittier method for donors where MACHINA migration graph may have some inconsistency with the data. 
        T = self.T
        cloneset = set()

        # query each node, recording reachable sites. 
        # missing: sites which variant WAS observed in, but UNREACHABLE by the clone. 
        # extra: sites which variant WAS NOT observed in, but REACHABLE by the clone. 
        candidates = []
        for node in sorted(list(nx.topological_sort(T))):
            clone = T.nodes[node]['clone']
            if clone not in self.valid:
                continue
            if clone in candidates:
                continue 
            reachable = set()
            reachable = reachable | set([T.nodes[n]['sample'] for n in T.nodes() if T.nodes[n]['clone']==clone]) # this clone's sample
            reachable = reachable | set(T.nodes[n]['sample'] for n in nx.descendants(T, node)) # descendent clones samples
            missing = observed - reachable
            extras = reachable - observed
            
            # Assume pidgeonhole. 
            # All observed samples must be reachable by clone. 
            # Extras may be simply due to low variant CCF in some samples.  
            # Especially true when the clone is in the primary prostate sample (and not truncal).
            if len(missing) == 0:
                candidates.append((clone, extras))

        # select candidate clones as those with minimal number of extra reachable samples. 
        min_extras = min([len(x[1]) for x in candidates])  
        for clone, extrasamples in candidates:
            if len(extrasamples) == min_extras:
                cloneset.add(clone)
        
        # Final choice: lowest candidate clone in the tree. 
        # This is because CCF gets lower as we descend the tree. 
        # Assumption is (as above), the issue was caused by low clone CCF in the primary prostate. 
        # This is why we don't 'observe' the variant in the prostate. 
        for node in sorted(list(nx.topological_sort(T)), reverse=True):
            clone = T.nodes[node]['clone']
            if clone in cloneset:
                return {clone}

        raise RuntimeError
    
    def _best_matching_clone_ccf(self, candidate_clones: list[str], candidate_ccfs: list[float], est_ccf: float) -> list:
        data = []
        for clone, c_ccf in zip(candidate_clones, candidate_ccfs):
            data.append((clone, c_ccf, est_ccf))
        df = pd.DataFrame.from_records(data=data, columns=['clone', 'clone_ccf', 'est_ccf'])
        df['diff'] = df['clone_ccf'] - df['est_ccf']
        df['diff'] = df['diff'].apply(lambda x: abs(x))
        df = df.sort_values('diff')
        return df.iloc[0].to_list()

class SNVassigner(CloneAssigner):

    def assign(self, coords: str, gene: str) -> Tuple:
        df = self.table[(self.table['coords']==coords) & (self.table['gene']==gene)].copy()
        info = self.assign_dpclust_direct(df, coords)
        if info is not None:
            return info
        return self.assign_site_parsimony_ccfs(df)
    
    def assign_dpclust_direct(self, df: pd.DataFrame, coords: str) -> Tuple|None:
        if coords not in self.asmts:
            return None
        
        clone = self.asmts[coords]
        sample = [self.T.nodes[n]['sample'] for n in self.T.nodes() if self.T.nodes[n]['clone']==clone][0]
        obs_ccfs = df[df['sample']==sample]['est_ccf'].values
        if len(obs_ccfs) == 0:
            return None 
        if len(obs_ccfs) >= 2:
            raise RuntimeError

        clone_ccf = self.ccfs.loc[clone, sample]
        obs_ccf = obs_ccfs[0]
        diff = abs(clone_ccf - obs_ccf)
        method = 'dpclust direct'
        return clone, clone_ccf, obs_ccf, diff, method
    
class INDELassigner(CloneAssigner):

    def assign(self, coords: str, gene: str) -> Tuple:
        df = self.table[(self.table['coords']==coords) & (self.table['gene']==gene)].copy()
        return self.assign_site_parsimony_ccfs(df)

class SVassigner(CloneAssigner):

    def assign(self, coords: str, gene: str) -> Tuple:
        df = self.table[(self.table['coords']==coords) & (self.table['gene']==gene)].copy()
        return self.assign_site_parsimony_earliest(df)

class CNAassigner(CloneAssigner):

    def assign(self, coords: str, gene: str) -> Tuple:
        df = self.table[(self.table['coords']==coords) & (self.table['gene']==gene)].copy()
        return self.assign_site_parsimony_ccfs(df)


def load_cmdline_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Merges all variant types into single file for a given PPCG donor.')
    parser.add_argument('--mutations', type=str, required=True, help='Path to merged variants (tsv).')
    parser.add_argument('--trees-dir', type=str, required=True, help='Path to directory containing machina tree data.')
    parser.add_argument('--clones-dir', type=str, required=True, help='Path to directory containing dpclust clone ccfs.')
    parser.add_argument('--outfile', type=str, required=True, help='Path to output file containing filtered variants.')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    main()
