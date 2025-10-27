
import argparse
import pandas as pd 
import numpy as np


def main() -> None:
    args = load_cmdline_args()
    
    # mutations
    muts = pd.read_csv(args.mutations, sep='\t', header=0)
    muts = muts[muts['vtype']!='CNApeak'].copy()

    # genesets
    gframe, sframe = load_genesets(args.genesets, args.selection)
    gframe = gframe[~gframe['gene'].str.startswith('MT-')].copy()
    sframe = sframe[~sframe['gene'].str.startswith('MT-')].copy()

    # genelists 
    background_genelist = gen_genelist(muts, gframe)

    # sizes 
    sizes = pd.read_csv(args.sizes, sep='\t', header=0)

    # mutation counts 
    mcounts = gen_mutation_counts(muts)

    # binary dtable
    matrix = gen_binary_matrix(muts, background_genelist)

    analyser = EnrichmentAnalyser(
        dtable=matrix,
        size_frame=sizes,
        gset_frame=sframe,
        mcounts_frame=mcounts,
        budget=args.budget
    )
    enrich = analyser.run()
    enrich = enrich.sort_values(['pval', 'factor'], ascending=[False, True])
    enrich.to_csv(args.outfile, sep='\t', index=False)


def load_genesets(genesets_path: str, selection_path: str):
    # load data 
    gframe = pd.read_csv(genesets_path, sep='\t', header=0)
    sel = pd.read_csv(selection_path, sep='\t', header=0)
    
    # make new genesets frame by subsetting selected
    all_genesets = sorted(list(gframe['geneset'].unique()))
    sel_genesets = sorted(list(sel[sel['valid']==True]['geneset'].unique()))
    assert len(sel_genesets) != 0
    assert len(set(sel_genesets) - set(all_genesets)) == 0
    sframe = gframe[gframe['geneset'].isin(sel_genesets)].copy()
    return gframe, sframe

def gen_genelist(muts: pd.DataFrame, gframe: pd.DataFrame):
    # background genelist
    mut_genes = set(muts['gene'].unique())
    full_gset_genes = set(gframe['gene'].unique())

    # remove mitochondrial genes
    mut_genes = set([x for x in mut_genes if not x.startswith('MT-')])
    full_gset_genes = set([x for x in full_gset_genes if not x.startswith('MT-')])
    
    # final genelist
    bkg_genes = set(mut_genes | full_gset_genes)

    print(f"mutated genes:    {len(mut_genes)}")
    print(f"geneset genes:    {len(full_gset_genes)}")
    print(f"union (genelist): {len(bkg_genes)}")
    return sorted(list(bkg_genes))

def gen_mutation_counts(muts: pd.DataFrame) -> pd.DataFrame:
    labels_l = ['seqvar', 'structvar', 'CNA']
    vclasses_l = [['SNV', 'INDEL'], ['SV'], ['CNA']]
    all_donors = sorted(list(muts['donor'].unique()))
    
    df = pd.DataFrame(index=all_donors)
    for label, vclasses in zip(labels_l, vclasses_l):
        mslice = muts[muts['vclass'].isin(vclasses)].copy()
        if label == 'seqvar':
            df[label] = mslice.groupby('donor')['coords'].nunique()
        else:
            df[label] = mslice.groupby('donor')['gene'].nunique()
    df = df.fillna(0).astype(int)
    return df

def gen_binary_matrix(muts: pd.DataFrame, genelist: list[str]) -> pd.DataFrame:
    all_donors = sorted(list(muts['donor'].unique()))
    df = pd.DataFrame(data=0, index=all_donors, columns=genelist)
    for rec in muts.itertuples():
        df.loc[rec.donor, rec.gene] = 1
    return df


class EnrichmentAnalyser:
    def __init__(
        self, 
        dtable: pd.DataFrame, 
        size_frame: pd.DataFrame, 
        gset_frame: pd.DataFrame,
        mcounts_frame: pd.DataFrame,
        budget: int=10
    ) -> None: 

        self.dtable  = dtable.T.copy()
        self.donors  = sorted(self.dtable.columns.to_list())
        self.genes   = sorted(self.dtable.index.to_list())
        self.mcounts = mcounts_frame
        self.budget  = budget
        print(self.dtable.iloc[:5, :5])
        
        self.seqvar_probs = self._get_gene_probabilities(self.genes, size_frame, field='cum_cds_len')
        self.structvar_probs = self._get_gene_probabilities(self.genes, size_frame, field='span')

        self.gset2genes = gset_frame.groupby('geneset')['gene'].agg(set).to_dict()
        self.gene2gsets = gset_frame.groupby('gene')['geneset'].agg(set).to_dict()
        for gene in self.genes:
            if gene not in self.gene2gsets:
                self.gene2gsets[gene] = set()
    
    def run(self) -> pd.DataFrame:
        obs_association = self.calc_observed_association()
        bkg_association = self.calc_background_association()
        results = self.calc_enrichment(obs_association, bkg_association)
        return results
    
    def calc_observed_association(self) -> pd.Series:
        return self._calc_hits(self.dtable)
    
    def calc_background_association(self) -> pd.DataFrame:
        assoc = pd.DataFrame()

        for i in range(self.budget):
            print(f'processed {i}/{self.budget} bootstrap samples...', end='\r')
            this_df = self._bootstrap_sample()
            this_assoc = self._calc_hits(this_df)
            assoc[i] = this_assoc
            if i % 10 == 0:
                assoc = assoc.copy()
        print(f'processed {i+1}/{self.budget} bootstrap samples... done.')
        assoc = assoc.fillna(0)
        assoc = assoc.astype(int)
        return assoc
    
    def calc_enrichment(self, obs_assoc: pd.Series, bkg_assoc: pd.DataFrame) -> pd.DataFrame:
        enrich_pvals = {}

        # pvals
        for k, v in obs_assoc.to_dict().items():
            pval = (bkg_assoc.loc[k] >= v).sum() / self.budget
            enrich_pvals[k] = min(pval, 1)

        # obs_assoc, mean background assoc
        enrich = pd.DataFrame.from_dict(enrich_pvals, orient='index')
        enrich.columns = ['pval']
        enrich['obs_donors'] = obs_assoc
        enrich['bkg_donors'] = bkg_assoc.mean(axis=1)
        enrich['factor'] = enrich['obs_donors'] / enrich['bkg_donors']
        enrich = enrich.reset_index()
        enrich = enrich.rename(columns={'index': 'geneset'})
        enrich = enrich.sort_values('pval')
        return enrich
    
    def _get_gene_probabilities(self, genelist: list[str], sframe: pd.DataFrame, field: str) -> list:
        sframe = sframe.set_index('gene')
        weights = pd.Series([sframe.loc[gene, field] for gene in genelist])
        probabilities = weights / weights.sum()
        return probabilities
    
    def _calc_hits(self, df: pd.DataFrame) -> pd.Series:
        data = []
        for gset, genes in self.gset2genes.items():
            ndonors = df.loc[list(genes)].sum().clip(upper=1).sum()
            data.append((gset, ndonors))
        df = pd.DataFrame.from_records(data=data, columns=['geneset', 'ndonors'])
        df = df.set_index('geneset')
        return df['ndonors']
    
    def _bootstrap_sample(self) -> pd.DataFrame:
        df = pd.DataFrame(data=0, index=self.genes, columns=self.donors)
        i = 0
        for donor in self.donors:
            n_seqvars = int(self.mcounts.loc[donor, 'seqvar'])  # type: ignore
            n_structvars = int(self.mcounts.loc[donor, 'structvar'])  # type: ignore
            n_cna = int(self.mcounts.loc[donor, 'CNA'])  # type: ignore

            genes_seqvar = np.random.choice(self.genes, size=n_seqvars, p=self.seqvar_probs, replace=True)
            genes_structvar = np.random.choice(self.genes, size=n_structvars, p=self.structvar_probs, replace=True)
            genes_cna = np.random.choice(self.genes, size=n_cna, replace=False)

            for genelist in [genes_seqvar, genes_structvar, genes_cna]:
                for gene in genelist:
                    df.loc[gene, donor] = 1

            i += 1
            if i % 20 == 0:
                df = df.copy()
        return df.copy()



def load_cmdline_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Performs differential mutation analysis between two cohorts.')
    parser.add_argument('--mutations', type=str, required=True, help='Path to mutations (tsv).')
    parser.add_argument('--genesets', type=str, required=True, help='Path to genesets (tsv).')
    parser.add_argument('--selection', type=str, required=True, help='Path to genesets selection (tsv).')
    parser.add_argument('--sizes', type=str, required=True, help='Path to gene sizes (tsv).')
    parser.add_argument('--budget', type=int, default=100, help='Number of bootstrap samples to use.')
    parser.add_argument('--outfile', type=str, required=True, help='Output file containing results (tsv).')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    main()