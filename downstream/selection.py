

import pandas as pd 

def select_genesets(gset_LUT: dict[str, list[str]], muts: pd.DataFrame):
    MAX_SIZE = 50
    MIN_DONORS = 10
    MIN_OVERLAP = 3
    MIN_COVERAGE = 0.2
    MAX_TOPGENE_RELIANCE = 0.5
    print()
    print(f"Total genesets in gmt:  {len(gset_LUT)}")

    # print('\nselecting genesets...')
    df = muts.copy()
    df = df[df['cohort']=='COMBI'].copy()
    # df = df[df['vclass']!='CNA'].copy()
    # df = muts.copy()

    data = []
    gcounts = df.groupby('gene')['donor'].nunique().sort_values(ascending=False)
    # print()
    # print(gcounts.head(10))

    for gset, genes in gset_LUT.items():
        genes = set(genes)
        dfslice = df[df['gene'].isin(genes)]
        shared = set(dfslice['gene'].unique()) & genes
        n_donors = dfslice['donor'].nunique()
        
        if len(shared) > 0:
            topgene = gcounts.loc[list(shared)].sort_values(ascending=False).index.to_list()[0]
        else:
            topgene = '__TOPGENE__'
        
        if n_donors > 0:
            n_donors_nontop = dfslice[dfslice['gene']!=topgene]['donor'].nunique()
            topgene_reliance = 1 - (n_donors_nontop / n_donors)
        else:
            n_donors_nontop = 0 
            topgene_reliance = 0

        data.append((
            gset,
            len(genes),
            len(shared),
            len(shared) / len(genes),
            n_donors,
            n_donors_nontop,
            topgene, 
            topgene_reliance
        ))
        
        # # geneset size 
        # if len(genes) > MAX_SIZE:
        #     ignore.add(gset)
        #     continue 
        
        # # overlap size
        # genes_gset = set(genes)
        # shared = genes_mut & genes_gset
        # if len(shared) < MIN_OVERLAP:
        #     ignore.add(gset)
        #     continue 

        # # overlap proportion (geneset coverage)
        # if len(shared) / len(genes_gset) < MIN_COVERAGE:
        #     ignore.add(gset)
        #     continue 

        # # donors affected
        # n_donors = df[df['gene'].isin(genes_gset)]['donor'].nunique()
        # if n_donors < MIN_DONORS:
        #     ignore.add(gset)

    df = pd.DataFrame.from_records(data, columns=['geneset', 'n_genes', 'n_genes_mut', 'coverage', 'n_donors', 'n_donors_nontop', 'topgene', 'topgene_reliance'])
    # print(df[df['topgene']=='TP53'].head(10))
    print()
    print(f"SIZE OK:         {(df['n_genes']<=MAX_SIZE).sum()/df.shape[0]*100:.1f}%")
    print(f"DONORS OK:       {(df['n_donors']>=MIN_DONORS).sum()/df.shape[0]*100:.1f}%")
    print(f"OVERLAP OK:      {(df['n_genes_mut']>=MIN_OVERLAP).sum()/df.shape[0]*100:.1f}%")
    print(f"COVERAGE OK:     {(df['coverage']>=MIN_COVERAGE).sum()/df.shape[0]*100:.1f}%")
    print(f"TOPGENE RELY OK: {(df['topgene_reliance']<=MAX_TOPGENE_RELIANCE).sum()/df.shape[0]*100:.1f}%")

    mask1 = df['n_genes']<=MAX_SIZE
    mask2 = df['n_donors']>=MIN_DONORS
    mask3 = df['n_genes_mut']>=MIN_OVERLAP
    mask4 = df['coverage']>=MIN_COVERAGE
    mask5 = df['topgene_reliance']<=MAX_TOPGENE_RELIANCE
    mask = mask1 & mask2 & mask3 & mask4 & mask5
    df = df[mask].copy()
    valid = set(df['geneset'].unique())
    gset_LUT = {k: v for k, v in gset_LUT.items() if k in valid}
    # print(f'{len(gset_LUT)} genesets retained.')
    print()
    print(f"Total genesets retained: {len(gset_LUT)}")
    return gset_LUT, df