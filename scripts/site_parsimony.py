import pandas as pd 
import networkx as nx 
from collections import defaultdict
from itertools import combinations

class SiteParsimonyAssigner:
    
    def __init__(self, T: nx.DiGraph, ccfs: pd.DataFrame) -> None:
        self._T = T
        self._ccfs = ccfs.drop('Cluster_Type', axis=1).clip(lower=0, upper=1).copy()
        self._sampat2clones_LUT = self._generate_sampat2clones_LUT()
        self._clone2ccfs_LUT = self._generate_clone2ccfs_LUT()

        # print('\n_sampat2clones_LUT')
        # for k, v in self._sampat2clones_LUT.items():
        #     print(k, v)
        # print('\n_clone2ccfs_LUT')
        # for k, v in self._clone2ccfs_LUT.items():
        #     print(k, v)

    def assign(self, obs_samples2ccfs: dict[str, float]) -> tuple[str, float, str]:
        assert len(obs_samples2ccfs) > 0
        obs_sams = sorted(list(obs_samples2ccfs.keys()))
        candidates = self._sampat2clones_LUT[tuple(obs_sams)]

        if len(candidates) == 0:
            raise RuntimeError
        
        # if len(candidates) == 1:
        #     if set(obs_samples2ccfs.keys()) == set(list()[])
        #     return list(candidates)[0], 'exact'
        
        scores = []
        for clone in candidates:
            dpc_samples2ccfs = self._clone2ccfs_LUT[clone]
            assert len(set(obs_samples2ccfs.keys()) - set(dpc_samples2ccfs.keys())) == 0
            score = 0
            statement = 'exact'
            for sample in dpc_samples2ccfs.keys():
                if sample not in obs_samples2ccfs:
                    score += 1 # penalty for variant not being witnessed in a sample it should 
                    statement = 'inexact'
                    continue
                score += abs(dpc_samples2ccfs[sample] - obs_samples2ccfs[sample])
            scores.append((clone, score, statement))
        
        scores = sorted(scores, key=lambda x: x[1])
        return scores[0]

    def _generate_clone2ccfs_LUT(self) -> dict:
        out = dict()
        for clone, row in self._ccfs.iterrows():
            out[clone] = row[row>0].to_dict()
        return out

    def _generate_sampat2clones_LUT(self) -> dict:
        out = dict()
        
        samples = sorted(self._ccfs.columns.to_list())
        itemsets = []
        for m in range(1, len(samples) + 1):
            itemsets.extend(combinations(samples, m))

        for comb in itemsets:
            out[comb] = self._get_clones(comb)
        return out

    def _get_clones(self, samples) -> set[str]:
        out = set()
        for node in nx.topological_sort(self._T):
            descendants = {node} | nx.descendants(self._T, node)
            reachable = set()
            for desc in descendants:
                reachable.update(self._T.nodes[desc]['samples'])
            if len(set(samples) - reachable) == 0:
                out.add(node)
        return out




if __name__ == '__main__':
    raise NotImplementedError


    # ccfs = load_dpclust(donor=PATIENT)
    # T = load_conipher(donor=PATIENT)
    # T = annotate_samples(T, ccfs)
    # samples = ccfs.drop('Cluster_Type', axis=1).columns.to_list()
    # draw_tree(T, ccfs)
    # print(ccfs)

    # assigner = SiteParsimonyAssigner(T, ccfs) 
    # tests = [
    #     {'PPCG0086a': 1.00, 'PPCG0086c': 1.00, 'PPCG0086d': 1.00, 'PPCG0086e': 1.00},
    #     {'PPCG0086a': 1.00, 'PPCG0086e': 1.00},
    #     {'PPCG0086c': 1.00, 'PPCG0086d': 1.00},
    #     {'PPCG0086d': 0.28},
    #     {'PPCG0086d': 0.95},
    #     {'PPCG0086e': 0.85},
    #     {'PPCG0086e': 0.11},
    #     {'PPCG0086a': 1.00, 'PPCG0086c': 1.00, 'PPCG0086d': 1.00},
    # ]
    # for test in tests:
    #     print()
    #     print(test)
    #     clone, score, is_exact = assigner.assign(obs_samples2ccfs=test)
    #     print(f"winner: {clone}, score={round(score, 2)}, {is_exact}")
