
import re
import sys
import numpy as np
import pandas as pd 
import networkx as nx 
from typing import Tuple

###############
### FILE IO ###
###############

def load_scna(filepath: str, allow_subclonal: bool) -> pd.DataFrame:
    df = pd.read_csv(filepath, sep='\t', header=0)
    df = df.rename(columns={'startpos': 'start', 'endpos': 'end'})
    df['chr'] = df['chr'].astype(str)
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)
    df = df[['chr', 'start', 'end', 'nMaj1_A', 'nMin1_A', 'frac1_A', 'nMaj2_A', 'nMin2_A', 'frac2_A']]
    
    data = []
    for idx, rec in df.iterrows():
        segments = []
        if not rec.isna()['nMaj2_A']:
            majority = max(rec.frac1_A, rec.frac2_A)
        else:
            majority = rec.frac1_A

        if rec.frac1_A > 0:
            segments.append((
                rec.chr, 
                rec.start, 
                rec.end, 
                rec.nMaj1_A, 
                rec.nMin1_A, 
                rec.nMaj1_A+rec.nMin1_A,
                rec.frac1_A / majority
            ))
        if not rec.isna()['frac2_A'] and rec.frac2_A > 0:
            segments.append((
                rec.chr, 
                rec.start, 
                rec.end, 
                rec.nMaj2_A, 
                rec.nMin2_A, 
                rec.nMaj2_A+rec.nMin2_A,
                rec.frac2_A / majority
            ))

        if len(segments)==2 and not allow_subclonal:
            continue 
        else:
            data += segments

    cnframe = pd.DataFrame.from_records(data, columns=['chr', 'start', 'end', 'nMaj', 'nMin', 'tcn', 'frac_ccf'])
    cnframe['chr'] = cnframe['chr'].astype('str')
    cnframe['start'] = cnframe['start'].astype('int')
    cnframe['end'] = cnframe['end'].astype('int')
    return cnframe


######################
### CCF ESTIMATION ###
######################

class CCFestimator:
    def __init__(self, purity: float, scna_path: str) -> None:
        self.purity = purity
        self.cnframe = load_scna(scna_path, allow_subclonal=True)

    def est_ccf(self, chrom: str, position: int, vaf: float) -> float:
        f = vaf 
        p = self.purity

        # get copy number information
        cnslice = self.cnframe[self.cnframe['chr']==chrom].copy()
        cnslice = cnslice[(cnslice['start']>=position) & (cnslice['end']<=position)].copy()
        if cnslice.shape[0] == 0:
            # assume normal copy number. 
            # tcn == 2 for autosomes, 1 for sex chromosomes
            segs = [[2.0, 1.0]] if chrom.isdigit() else [[1.0, 1.0]]
        else:
            segs = [[float(rec.tcn), rec.frac_ccf] for rec in cnslice.itertuples()] # type: ignore

        data = []
        for nt, frac_ccf in segs:
            m = max(1, round(( f/p ) * ( p*nt + 2*(1-p) )))
            ccf = ( f/(m*p) ) * ( p*nt + 2*(1-p) )
            diff = abs(frac_ccf-ccf)
            data.append((round(ccf, 2), diff))
        data.sort(key=lambda x: x[-1])
        return min(1, data[0][-2])
    




#################
### FILTERING ###
#################




# def load_dpclust_ccfs(filepath: str) -> pd.DataFrame:
#     df = pd.read_csv(filepath, header=0)
#     df = df.rename(columns={'Cluster': 'clone'})
#     df.columns = [x.replace('_DNA', '') for x in df.columns]
#     df['clone'] = df['clone'].astype(str)
#     df['clone'] = df['clone'].apply(lambda x: x.replace('_', ''))
#     df = df.set_index('clone')
#     return df

# def load_dpclust_assignments(filepath: str) -> pd.Series:
#     df = pd.read_csv(filepath, header=0)
#     df.columns = [x.lower() if not x.startswith('PPCG') else x for x in df.columns]
#     df.columns = [x.replace('_DNA', '') for x in df.columns]
#     df['cluster'] = df['cluster'].astype(str)
#     df['cluster'] = df['cluster'].apply(lambda x: x.replace('_', ''))
#     df['chr'] = df['chr'].astype(str)
#     df['pos'] = df['pos'].astype(str)
#     df['coords'] = df['chr'] + ':' + df['pos']
#     df = df.set_index('coords')
#     df = df['cluster']
#     return df

# def load_clonetree(tree_labels: str, tree_edges: str, tree_dotfile: str, sam2site: str) -> nx.DiGraph:
#     PATTERN_DOTPROP = r'(\d+) \[.*label="([\d\.]+%[^"]*)"\]'
#     PATTERN_DOTEDGE = r'(\d+|X) -> (\d+) \[.*label="([\w^]+)"\]'
#     PATTERN_CLONE = r'(\d+).*'

#     # node sites
#     with open(tree_labels, 'r') as fp:
#         lines = fp.readlines()
#         lines = [ln.strip().split(' ') for ln in lines]
#         site_map: dict[str, str] = {k: v for k, v in lines}

#     # leaf node proportions
#     with open(tree_dotfile, 'r') as fp:
#         lines = fp.readlines()
#         lines = [ln.strip() for ln in lines]
        
#         # get leaf node proportions
#         mprops_map = {}
#         for ln in lines:
#             m = re.match(PATTERN_DOTPROP, ln)
#             if m is not None:
#                 mnode, props_str = m.group(1), m.group(2)
#                 props_l = [float(x.strip('%')) for x in props_str.split('\\n')]
#                 maxprop = max(props_l)
#                 mprops_map[mnode] = maxprop
        
#         cprops_map = {}
#         for ln in lines:
#             m = re.match(PATTERN_DOTEDGE, ln)
#             if m is not None:
#                 mnode, node = m.group(2), m.group(3)
#                 if mnode in mprops_map:
#                     cprops_map[node] = mprops_map[mnode]

#     # graph structure from edgelist
#     smapper = pd.read_csv(sam2site, sep='\t', header=0)
#     site2sample = smapper.set_index('site')['sample'].to_dict()
#     T = nx.DiGraph()
#     with open(tree_edges, 'r') as fp:
#         line = fp.readline().strip()
#         while line:
#             src_label, dest_label = line.split(' ')
#             src_m = re.match(PATTERN_CLONE, src_label)
#             dest_m = re.match(PATTERN_CLONE, dest_label)
#             assert src_m
#             assert dest_m
#             src_clone = src_m.group(1)
#             src_site = site_map[src_label]
#             src_sample = site2sample[src_site]

#             dest_clone = dest_m.group(1)
#             dest_site = site_map[dest_label]
#             dest_sample = site2sample[dest_site]

#             src_prop = None if src_label not in cprops_map else cprops_map[src_label]
#             dest_prop = None if dest_label not in cprops_map else cprops_map[dest_label]
#             T.add_node(src_label, clone=src_clone, site=src_site, sample=src_sample, prop=src_prop)
#             T.add_node(dest_label, clone=dest_clone, site=dest_site, sample=dest_sample, prop=dest_prop)
#             T.add_edge(src_label, dest_label)

#             line = fp.readline().strip()

#     for node in T.nodes():
#         if T.out_degree(node) == 0:
#             assert T.nodes[node]['prop'] is not None

#     return T

# def modify_structure(T: nx.DiGraph) -> nx.DiGraph:
#     # remove low CCF trunk metastatic events (likely error, they don't make sense, eg PPCG0180)

#     T = _remove_lowccf_trunk_metevents(T)

#     # make node sites labelling consistent / sensible. MACHINA has odd output in some cases. 
#     for node in list(T.nodes()):
#         parents = list(T.predecessors(node))
#         assert len(parents) <= 1
#         if _should_reassign_parent(T, node, parents):
#             T = _do_reassign_parent(T, node, parents[0])
    
#     # why the fuck is this here again? doing it 2 times? wtf? 
#     for node in list(T.nodes()):
#         parents = list(T.predecessors(node))
#         assert len(parents) <= 1
#         if _should_reassign_parent(T, node, parents):
#             print(f'WARNING: PARENT REASSIGNMENT IN REPEATED LOOP')
#             T = _do_reassign_parent(T, node, parents[0])
    
#     # standardises format for metastatic events.
#     for node in list(T.nodes()):
#         parents = list(T.predecessors(node))
#         assert len(parents) <= 1
#         if _should_add_placeholder(T, node, parents):
#             T = _do_add_placeholder(T, node, parents[0])

#     T = _propagate_proportions(T)
#     T = _prune_leaves(T)

#     # validation
#     metevents = _get_seedclones(T)
#     if len(metevents) == 0:
#         raise RuntimeError
    
#     return T

# def _propagate_proportions(T: nx.DiGraph) -> nx.DiGraph:
#     for parent in list(nx.topological_sort(T))[::-1]:
#         if T.out_degree(parent) > 0: # internal and root
#             prop = 0
#             for child in T.successors(parent):
#                 if T.nodes[parent]['site'] == T.nodes[child]['site']:
#                     prop += T.nodes[child]['prop']
#             assert prop <= 102
#             prop = min(100, prop)
#             T.nodes[parent]['prop'] = prop 
    
#     return T

# def _remove_lowccf_trunk_metevents(T: nx.DiGraph) -> nx.DiGraph:
#     THRESHOLD = 5   # 5% estimated sample CCF
#     root = list(nx.topological_sort(T))[0]
#     leaves = [n for n in T.nodes() if len(list(T.successors(n)))==0]
#     for leaf in leaves:
#         prop = T.nodes[leaf]['prop']
#         if prop is None:
#             print(leaf)
#     leaves = [n for n in leaves if T.nodes[n]['prop']<THRESHOLD]
#     # print('Removing low CCF leaves ---')
#     for leaf in leaves:
#         parents_l = list(T.predecessors(leaf))
#         assert len(parents_l) <= 1
#         parent = parents_l[0]
#         if parent == root and T.nodes[leaf]['site'] != 'prostate':
#             # print(f"Removing {leaf}\t({T.nodes[leaf]['prop']}%)")
#             T.remove_node(leaf)
#     return T

# def _get_seedclones(T: nx.DiGraph) -> list[Tuple]:
#     events = list()
#     for n1, n2 in T.edges():
#         if T.nodes[n1]['clone'] != T.nodes[n2]['clone']:
#             continue
#         if T.nodes[n1]['site'] == T.nodes[n2]['site']:
#             continue
#         events.append((n1, n2))
#     return events

# def _should_reassign_parent(T: nx.DiGraph, node: str, parents: list[str]) -> bool:
#     # checking there is a parent & grandparent
#     if len(parents) != 1:
#         return False
#     parent = parents[0]
#     grandparents = list(T.predecessors(parent)) 
#     if len(grandparents) != 1:
#         return False
#     grandparent = grandparents[0]

#     # checking the parent only has 1 child
#     peers = list(T.successors(parent))
#     assert len(peers) >= 1
#     if len(peers) != 1:
#         return False
    
#     # checking the grandparent and parent have the same site, 
#     # and the child node has a different site
#     grand_site = T.nodes[grandparent]['site']
#     parent_site = T.nodes[parent]['site']
#     child_site = T.nodes[node]['site']
#     if grand_site != parent_site:
#         return False
#     if parent_site == child_site:
#         return False

#     return True

# def _do_reassign_parent(T: nx.DiGraph, node: str, parent: str) -> nx.DiGraph:
#     T.nodes[parent]['site'] = T.nodes[node]['site']
#     T.nodes[parent]['sample'] = T.nodes[node]['sample']
#     return T

# def _should_add_placeholder(T: nx.DiGraph, node: str, parents: list[str]) -> bool:
#     PATTERN = r'^(\d+)\^([a-zA-Z]\w+)$'
#     if len(parents) != 1:
#         return False
#     parent = parents[0]
#     # print(f"parent: {parent}, {T.nodes[parent]['site']}")
#     # print(f"child: {node}, {T.nodes[node]['site']}")
#     # print()
#     if T.nodes[node]['site'] == T.nodes[parent]['site']:
#         return False
#     if re.match(PATTERN, node) is not None:
#         return False
#     return True

# def _do_add_placeholder(T: nx.DiGraph, node: str, parent: str) -> nx.DiGraph:
#     PATTERN = r'^(\d+)_([a-zA-Z]\w+)$'
#     site = T.nodes[node]['site']
#     sample = T.nodes[node]['sample']
#     if re.match(PATTERN, node) is not None:
#         new_node = node.replace('_', '^', 1)
#         new_clone = re.match(PATTERN, node).group(1)  # type: ignore
#     else:
#         new_node = f"{parent}^{site}"
#         new_clone = parent
#     # print(f'adding placeholder for node {node}: {new_clone}')
#     T.add_node(new_node, clone=new_clone, site=site, sample=sample)
#     T.add_edge(parent, new_node)
#     T.add_edge(new_node, node)
#     T.remove_edge(parent, node)
#     return T

# def _prune_leaves(T: nx.DiGraph) -> nx.DiGraph:
#     # pruning leaves
#     leaves = [n for n in T.nodes() if T.out_degree(n) == 0]
#     for leaf in leaves:
#         # print(f'removing leaf: {leaf}')
#         T.remove_node(leaf)
#     return T


# ########################
# ### CLONE ASSIGNMENT ###
# ########################

# class SiteParsimonyAssigner:
#     def __init__(self, donor: str) -> None:

      
#     def assign_cna(self, seginfo: list[Tuple[str, list[SCNSegment]]]) -> Tuple:
#         """
#         The seginfo variable holds 1+ SCNA segment, each with 1+ SNAallele per sample. 
#         Reason for multiple segments: 
#           - A given sample can contain multiple different SCNA segments in the query interval. 
#           - eg. query interval chr1:1-10_000, segment1=chr1:1-2_000, segment2=chr1:2_001-10_000 
#         Reason for multiple alleles:
#           - Segments are subclonal. 
#           - Truncal segments will have only 1 allele a1. 
#           - Subclonal segments will have 2 alleles: a1 & a2. 
#         """
#         obs_samples = sorted([x[0] for x in seginfo])
#         candidate_sample, candidate_clones, candidate_ccfs, method = self._parsimonious_clones(obs_samples)
        
#         # only 1 parsimonious clone. Simple. 
#         if len(candidate_clones) == 1:
#             return candidate_clones.pop(), None, None, None
        
#         # 2+ parsimonious clones. Use CCF similarity. 
#         # First, get the segments which belong to the parsimonious sample. 
#         selected = [(sam, seg) for sam, seg in seginfo if sam==candidate_sample]
#         assert len(selected) == 1
#         segments = selected[0][1]

#         # Due to the comments in this function's docstring, we will use the truncal segment (if exists). 
#         for seg in segments:
#             if seg.a1.ccf == 1 and seg.a2 is None:
#                 # this function is somewhat unnecessary but avoids more code. 
#                 return self._best_matching_clone_ccf(candidate_clones, candidate_ccfs, [1.0])

#         # If no truncal segment, we consider *each* SCNA allele for *each* segment. 
#         # The candidate clone which has most similar CCF to one of these SCNA allele CCFs is selected. 
#         ccfs = []
#         for seg in segments:
#             assert seg.a2 is not None
#             ccfs.append(seg.a1.ccf)
#             ccfs.append(seg.a2.ccf)
#         return self._best_matching_clone_ccf(candidate_clones, candidate_ccfs, ccfs)

#     def assign_general(self, obs_samples: list[str], obs_vafs: list[float], chrom: str, position: int, svclass: str|None=None) -> Tuple:
#         sample, candidate_clones, candidate_ccfs, method = self._parsimonious_clones(obs_samples)

#         # only 1 parsimonious clone. Simple. 
#         if len(candidate_clones) == 1:
#             return candidate_clones.pop(), None, None, None
        
#         # 2+ parsimonious clones. Use CCF similarity. 
#         # If this is a SNV/INDEL/SV mutation, the obs_vafs variable holds 1 VAF per sample.
#         # Now we know *which* sample the variant arose in, we select the VAF for this sample. 
#         # This VAF is converted into a CCF before comparing against candidate clones.
#         obs_vaf = obs_vafs[obs_samples.index(sample)]
#         sub_ccfs = self._calc_ccfs_subclonal(sample, obs_vaf, chrom, position, svclass)
#         return self._best_matching_clone_ccf(candidate_clones, candidate_ccfs, sub_ccfs)
    
#     def _parsimonious_clones(self, obs_samples: list[str]) -> Tuple:
#         """
#         Assigns mutation to best fitting clone using observed pattern of anatomical sites. 
#         """
#         obssites = set([self.sam2site_lut[s] for s in obs_samples])
#         cloneset = self._parsimonious_clones_exact(obssites)
#         if len(cloneset) >= 1:
#             c_clones = cloneset
#             method = 'EXACT'
#         else:
#             c_clones = self._parsimonious_clones_inexact(obssites)
#             method = 'INEXACT'
#         c_sites = set([self.T.nodes[n]['site'] for n in c_clones])
        
#         # ensure 1+ candidate clones
#         assert len(c_clones) >= 1
        
#         # ensure single site/sample for parsimonious clones 
#         if len(c_sites) != 1:
#             raise RuntimeError('2+ anatomical sites.')

#         sample = self.sam2site_lut_r[c_sites.pop()]
#         candidate_clones = sorted(list(c_clones))
#         candidate_ccfs = [float(self.dpclust_ccfs.loc[c, sample]) for c in candidate_clones]
#         return sample, candidate_clones, candidate_ccfs, method 

#     def _parsimonious_clones_exact(self, obssites: set[str]) -> set[str]:
#         # correct method
#         T = self.T
#         cloneset = set()
#         for node in nx.topological_sort(T):
#             clone = T.nodes[node]['clone']
#             if clone not in self.valid:
#                 continue
#             # this site
#             desc_sites = set([T.nodes[n]['site'] for n in T.nodes() if T.nodes[n]['clone']==clone])
#             # reachable sites
#             desc_sites = desc_sites | set(T.nodes[n]['site'] for n in nx.descendants(T, node))
#             if desc_sites == obssites:
#                 cloneset.add(clone)
#         return cloneset

#     def _parsimonious_clones_inexact(self, obssites: set[str]) -> set[str]:
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
#             desc_sites = set([T.nodes[n]['site'] for n in T.nodes() if T.nodes[n]['clone']==clone])
#             desc_sites = desc_sites | set(T.nodes[n]['site'] for n in nx.descendants(T, node))
#             missing = obssites - desc_sites
#             extras = desc_sites - obssites
            
#             # Assume pidgeonhole. 
#             # All observed samples must be reachable by clone. 
#             # Extras may be simply due to low variant CCF in some samples.  
#             # Especially true when the clone is in the primary prostate sample (and not truncal).
#             if len(missing) == 0:
#                 candidates.append((clone, extras))

#         # select candidate clones as those with minimal number of extra reachable sites. 
#         min_extras = min([len(x[1]) for x in candidates])  
#         for clone, extrasites in candidates:
#             if len(extrasites) == min_extras:
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
        
#     def _calc_ccfs_subclonal(
#         self, 
#         sample: str, 
#         vaf: float, 
#         chrom: str, 
#         position: int, 
#         svclass: str|None
#         ) -> list[float]: 
#         # rename to provide clarity with formula 
#         f = vaf
#         p = self.sam2purity_lut[sample]
#         nt_l = self._local_cns_at(sample, chrom, position)
        
#         # calculate ccf for each subclonal segment at this location
#         ccfs = []
#         for nt in nt_l:
#             # calculate multiplicity
#             # m = ( f/p ) * ( p*nt + 2*(1-p) )
#             if svclass == 'DUP':
#                 m = 3
#             else:
#                 m = 2

#             # calculate ccf
#             ccf = ( f/(m*p) ) * ( p*nt + 2*(1-p) )
#             ccfs.append(ccf)
#             # print(f"f={f}, p={p}, nt={nt}, m={m:.2f}, ccf={ccf:.2f}")
        
#         return ccfs
        
#     def _local_cns_at(self, sample: str, chrom: str, position: int) -> list[int]:
#         segments = self.sam2segs_lut[sample][chrom]
#         for seg in segments: 
#             if seg.start <= position and seg.end >= position:
#                 return [seg.a1.tot_cn] if seg.a2 is None else  [seg.a1.tot_cn, seg.a2.tot_cn]
#         # TODO: No segment, assume normal? wtf? what about WGD?
#         return [2]
#         # if self.wgdmapper.sample(sample):
#         #     return [4]
#         # return [2]
    
#     def _best_matching_clone_ccf(self, candidate_clones: list[str], candidate_ccfs: list[float], est_ccfs: list[float]) -> Tuple:
#         data = []
#         for clone, c_ccf in zip(candidate_clones, candidate_ccfs):
#             for e_ccf in est_ccfs:
#                 data.append((clone, c_ccf, e_ccf))
#         df = pd.DataFrame.from_records(data=data, columns=['clone', 'clone_ccf', 'est_ccf'])
#         df['diff'] = df['clone_ccf'] - df['est_ccf']
#         df['diff'] = df['diff'].apply(lambda x: abs(x))
#         df = df.sort_values('diff')
#         # print(df)
#         clone, clone_ccf, est_ccf, diff = df.iloc[0].to_list()
#         return clone, clone_ccf, est_ccf, diff
        

