
import pandas as pd 
import networkx as nx
import matplotlib.pyplot as plt
from itertools import combinations
from networkx.drawing.nx_agraph import graphviz_layout

######################
### Site-Parsimony ###
######################

class SiteParsimonyAssigner:
    
    def __init__(self, ccfs_path: str, tree_path: str) -> None:
        self.ccfs = self._load_dpclust_ccfs(ccfs_path)
        self.tree = self._load_conipher_tree(tree_path)
        self._sampat2clones_LUT = self._generate_sampat2clones_LUT()
        self._clone2ccfs_LUT = self._generate_clone2ccfs_LUT()

        # print('\n_sampat2clones_LUT')
        # for k, v in self._sampat2clones_LUT.items():
        #     print(k, v)
        # print('\n_clone2ccfs_LUT')
        # for k, v in self._clone2ccfs_LUT.items():
        #     print(k, v)
    
    def assign(self, obs_samples2ccfs: dict[str, float]) -> tuple[str, float, bool]:
        assert len(obs_samples2ccfs) > 0
        obs_sams = sorted(list(obs_samples2ccfs.keys()))
        candidates = self._sampat2clones_LUT[tuple(obs_sams)]

        if len(candidates) == 0:
            raise RuntimeError
        
        scores = []
        for clone in candidates:
            dpc_samples2ccfs = self._clone2ccfs_LUT[clone]
            assert len(set(obs_samples2ccfs.keys()) - set(dpc_samples2ccfs.keys())) == 0
            score = 0
            exact = True
            for sample in dpc_samples2ccfs.keys():
                if sample not in obs_samples2ccfs:
                    score += 1 # penalty for variant not being witnessed in a sample it should 
                    exact = False
                    continue
                score += abs(dpc_samples2ccfs[sample] - obs_samples2ccfs[sample])
            scores.append((clone, score, exact))
        
        scores = sorted(scores, key=lambda x: x[1])
        return scores[0]
    
    def _load_dpclust_ccfs(self, filepath: str, minccf: float=0.05) -> pd.DataFrame:
        df = pd.read_csv(filepath, sep=',', header=0)
        df['Cluster'] = df['Cluster'].apply(lambda x: str(x).replace('_', '').split('.')[0])
        df = df.set_index('Cluster')
        df.columns = [x.replace('_DNA', '') for x in df.columns]
        df = df[[x for x in df.columns if x.startswith('PPCG')]].copy()
        for col in df.columns:
            mask = df[col]<minccf
            df.loc[mask, col] = 0.00
        df = df.clip(lower=0, upper=1)
        return df.copy()
    
    def _load_conipher_tree(self, filepath: str) -> nx.DiGraph:
        patient = self.ccfs.columns.to_list()[0][:8]
        T = self._load_tree1(filepath)
        
        # hacky hot-fix for PPCG0435 where clone 11 was merged into clone 4 (trunk)
        if patient == 'PPCG0435':
            T.remove_node('11')
            
        T = self._annotate_samples(T)
        return T

    def _annotate_samples(self, T: nx.DiGraph) -> nx.DiGraph:
        df = self.ccfs.T.copy()
        for node in T.nodes():
            present = df[df[node]>0].index.to_list()
            T.nodes[node]['samples'] = set(present)
        return T

    def _load_tree1(self, filepath: str) -> nx.DiGraph:

        trees = {}
        
        with open(filepath, 'r') as fp:
            line = fp.readline().strip()
            line = fp.readline().strip()
            assert line == '# tree 1'

            name = line.strip('# ')
            T = nx.DiGraph()
            line = fp.readline().strip()
            while line:
                if line.startswith('#'):
                    trees[name] = T
                    name = line.strip('# ')
                    T = nx.DiGraph()
                else:
                    parent, child = line.split('\t')
                    parent = parent.replace('_', '')
                    child = child.replace('_', '')
                    T.add_edge(parent, child)
                line = fp.readline().strip()
            trees[name] = T
        
        return trees['tree 1']

    def _generate_clone2ccfs_LUT(self) -> dict:
        out = dict()
        for clone, row in self.ccfs.iterrows():
            out[clone] = row[row>0].to_dict()
        return out

    def _generate_sampat2clones_LUT(self) -> dict:
        out = dict()
        
        samples = sorted(self.ccfs.columns.to_list())
        itemsets = []
        for m in range(1, len(samples) + 1):
            itemsets.extend(combinations(samples, m))

        for comb in itemsets:
            out[comb] = self._get_clones(comb)
        return out

    def _get_clones(self, samples) -> set[str]:
        out = set()
        for node in nx.topological_sort(self.tree):
            descendants = {node} | nx.descendants(self.tree, node)
            reachable = set()
            for desc in descendants:
                reachable.update(self.tree.nodes[desc]['samples'])
            if len(set(samples) - reachable) == 0:
                out.add(node)
        return out
    
    def draw_tree(self) -> None:
        T = self.tree
        df = self.ccfs.copy()
        samples = df.columns.to_list()
        
        i = 0
        ncols = 4 if len(samples) >= 4 else len(samples) % 4
        nrows = len(samples) // 4 
        if len(samples) % 4 != 0:
            nrows += 1
        
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(4*ncols, 4*nrows))
        for sample in samples:
            if nrows > 1:
                row = i // 4
                col = i % 4 
                ax = axes[row, col]
            else:
                col = i % 4 
                ax = axes[col]

            pos = graphviz_layout(T, prog='dot')
            present = df[df[sample]>0].index.to_list()
            nl_present = [n for n in T.nodes() if n in present]
            nl_missing = [n for n in T.nodes() if n not in present]
            nx.draw_networkx_nodes(T, pos, nodelist=nl_present, node_color='skyblue', node_size=500, ax=ax)
            nx.draw_networkx_nodes(T, pos, nodelist=nl_missing, node_color='salmon', node_size=500, ax=ax)
            nx.draw_networkx_labels(T, pos, ax=ax)
            nx.draw_networkx_edges(T, pos, ax=ax)
            ax.set_title(sample)
            i += 1
        plt.show()
        plt.close()


if __name__ == '__main__':
    
    DONOR = 'PPCG0086'
    INDIR_DPCLUST = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/angel/26.03.2026/ccf/Clonal_trees'
    INDIR_CONIPHER = '/home/grace/work/PPCG_DifferentialGenesetMutation/data/angel/06.03.2026/conipher_trees'
    
    tree_path = f"{INDIR_CONIPHER}/{DONOR}_conipher_tree/allTrees.txt"
    ccfs_path = f"{INDIR_DPCLUST}/{DONOR}_Cluster_CCFs.csv"
    assigner = SiteParsimonyAssigner(ccfs_path=ccfs_path, tree_path=tree_path)
    # assigner.draw_tree()
    print(assigner.ccfs)
    tests = [
        {'PPCG0086a': 1.00, 'PPCG0086c': 1.00, 'PPCG0086d': 1.00, 'PPCG0086e': 1.00},
        {'PPCG0086a': 1.00, 'PPCG0086e': 1.00},
        {'PPCG0086c': 1.00, 'PPCG0086d': 1.00},
        {'PPCG0086a': 1.00, 'PPCG0086c': 1.00},
        {'PPCG0086d': 0.28},
        {'PPCG0086d': 0.95},
        {'PPCG0086e': 0.85},
        {'PPCG0086e': 0.11},
        {'PPCG0086a': 1.00, 'PPCG0086c': 1.00, 'PPCG0086d': 1.00},
    ]
    for test in tests:
        clone, score, is_exact = assigner.assign(obs_samples2ccfs=test)
        meth = 'site-parsimony (exact)' if is_exact else 'site-parsimony (inexact)'
        print(f"winner={clone}, score={round(score, 2)}, method={meth}, obs={test}")

