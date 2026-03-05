
import argparse
import pandas as pd 
import matplotlib.pyplot as plt 
import PyComplexHeatmap as pycomp
from scipy import stats 

pd.options.display.float_format = '{:.2f}'.format

def main() -> None:
    args = load_cmdline_args()

    ### LOAD DATA 
    pmuts = pd.read_csv(args.posmuts, sep='\t', header=0)
    nmuts = pd.read_csv(args.negmuts, sep='\t', header=0)
    gframe = pd.read_csv(args.genesets, sep='\t', header=0)
    sel = pd.read_csv(args.selection, sep='\t', header=0)
    diff = pd.read_csv(args.diffmut, sep='\t', header=0)
    prich = pd.read_csv(args.posenrich, sep='\t', header=0)
    nrich = pd.read_csv(args.negenrich, sep='\t', header=0)
    sel_genesets = sorted(list(sel[sel['valid']==True]['geneset'].unique()))
    diff = diff.set_index('geneset')
    prich = prich.set_index('geneset')
    nrich = nrich.set_index('geneset')
    # prich['pval_fdr'] = stats.false_discovery_control(prich['pval'].to_list())
    # nrich['pval_fdr'] = stats.false_discovery_control(nrich['pval'].to_list())

    print(gframe['geneset'].nunique())
    print(sel['geneset'].nunique())
    print(len(sel_genesets))

    n_donors_pos = pmuts['donor'].nunique()
    n_donors_neg = nmuts['donor'].nunique()
    print(f"pos donors: {n_donors_pos}")
    print(f"neg donors: {n_donors_neg}")

    df = pd.DataFrame(index=sel_genesets)
    df['n_genes'] = gframe[gframe['geneset'].isin(sel_genesets)].groupby('geneset')['gene'].nunique()
    df['obs_treat'] = diff['posclass donors']
    df['obs_control'] = diff['negclass donors']
    df['bkg_treat'] = prich['background']
    df['bkg_control'] = nrich['background']
    df['fisher_pval'] = diff['fisher_pval']
    df['logit_pval'] = diff['logit_pval']
    # df['diffmut_pval'] = diff['logit_pval']
    df['bkg_pval_treat'] = prich['pval']
    df['bkg_pval_control'] = nrich['pval']
    df['diffmut_factor'] = df['obs_treat'] / df['obs_control']
    df['bkg_factor_treat'] = prich['factor']
    df['bkg_factor_control'] = nrich['factor']
    
    # # sorting
    # mask = df['obs_control']!=0
    # df.loc[mask, 'sort_factor'] = df.loc[mask, 'diffmut_factor'].apply(lambda x: round(x))
    # df.loc[~mask, 'sort_factor'] = df['obs_treat']
    # df['sort_logit'] = df['diffmut_pval'].apply(lambda x: round(x, 2))
    # df = df.sort_values(['sort_logit', 'sort_factor'], ascending=[True, False])
    # df = df.drop(['sort_logit', 'sort_factor'], axis=1)

    # filter
    def is_significant(row: pd.Series) -> bool:
        # # has to be significantly enriched in positive cohort. 
        if row['fisher_pval_FDR'] > 0.1:
            return False 
        return True
        # if row['bkg_pval_treat'] > 0.1:
        #     return False 
        # if row['bkg_pval_control'] <= 0.1:
        #     return False 

        # # has to be enriched above background in positive cohort. 
        # if row['bkg_factor_treat'] < 1:
        #     return False 
        # if row['bkg_pval_treat'] > 0.05:
        #     return False
        # # can't be observed in negative cohort significantly below expected background rate
        # if row['bkg_factor_control'] < 1 and row['bkg_factor_control'] <= 0.05:
        #     return False 
        # return True 

    # filtering
    df['fisher_pval_FDR'] = stats.false_discovery_control(df['fisher_pval'].to_list())
    df['sig'] = df.apply(is_significant, axis=1)
    outfile = f"{args.posclass}.full.tsv"
    df = df.reset_index()
    df = df.sort_values('fisher_pval')
    df = df[[x for x in df.columns if x!='index']+['index']]
    df.to_csv(outfile, sep='\t', index=False, float_format='%.3f')
    return

    # df_full = df.copy()
    # df_sig = df.copy()
    # df_sig['valid'] = df_sig.apply(is_valid, axis=1)
    # df_sig = df_sig[df_sig['valid']==True].copy()
    # df_sig = df_sig.drop('valid', axis=1)
    # df_sig['diffmut_pval_fdr'] = stats.false_discovery_control(df_sig['diffmut_pval'].values)
    # df_sig = df_sig[df_sig['diffmut_pval_fdr']<=0.1].copy()
    # df_sig = df_sig.drop('diffmut_pval', axis=1)

    # # write to file
    # rframe = df_sig.copy()
    # outfile_full = f"{args.posclass}.full.tsv"
    # outfile_sig = f"{args.posclass}.sig.tsv"
    # df_full = df_full.reset_index()
    # df_sig = df_sig.reset_index()
    # df_full = df_full[[x for x in df_full.columns if x!='index']+['index']]
    # df_sig = df_sig[[x for x in df_sig.columns if x!='index']+['index']]
    # df_full.to_csv(outfile_full, sep='\t', index=False, float_format='%.3f')
    # df_sig.to_csv(outfile_sig, sep='\t', index=False, float_format='%.3f')
    
    # merge mutations
    POSCLASS = f"{args.posclass} ({n_donors_pos} donors)"
    NEGCLASS = f"Others ({n_donors_neg} donors)"
    pmuts['cohort'] = POSCLASS
    nmuts['cohort'] = NEGCLASS
    muts = pd.concat([pmuts, nmuts], ignore_index=True)
    print(muts['cohort'].value_counts())
    print(muts.head())

    # for plotting later later. 
    CMAPPER = {
        'SV': '#9467bd',
        'SNV': '#1f77b4',
        'INDEL': '#ff7f0e',
        'CNA↑': '#2ca02c',
        'CNA↓': '#d62728',
    }
    mask_up = muts['vtype'].isin(['CNA↑'])
    mask_dn = muts['vtype'].isin(['CNA↓', 'CNA↓↓'])
    muts.loc[mask_up, 'vclass'] = 'CNA↑'
    muts.loc[mask_dn, 'vclass'] = 'CNA↓'
    print()
    print(muts['vclass'].value_counts())

    # Columns are donors mapped to cohorts. 
    col_split = muts.drop_duplicates(subset=['donor'])[['donor', 'cohort']]
    col_split = col_split.set_index('donor')
    print(col_split.head())
    print(col_split.tail())

    # TMB
    tmb_frame = muts.groupby('donor')['vclass'].value_counts().unstack().fillna(0).astype(int)
    tmb_frame['total'] = tmb_frame.sum(axis=1)
    tmb_frame = tmb_frame.sort_values('total', ascending=False)
    tmb_frame['cohort'] = col_split['cohort']
    print()
    print(tmb_frame.head(10))
    tmb_frame = tmb_frame.drop(['total', 'cohort'], axis=1)

    # OncoPrints
    top_n = rframe.head(args.max_plots).index.to_list()
    for geneset in top_n:
        outfile = f"oncoprints/{geneset}.png"
        unstacked = gen_unstacked(muts, gframe, geneset=geneset, posclass=POSCLASS, negclass=NEGCLASS, cmapper=CMAPPER)
        render_oncoprint(unstacked, col_split, tmb_frame, posclass=POSCLASS, negclass=NEGCLASS, cmapper=CMAPPER, filepath=outfile)

def gen_unstacked(mtable: pd.DataFrame, gframe: pd.DataFrame, geneset: str, posclass: str, negclass: str, cmapper: dict) -> pd.DataFrame:
    idents = []
    membership_l = [posclass, negclass]
    donors_pos = sorted(list(mtable[mtable['cohort']==posclass]['donor'].unique()))
    donors_neg = sorted(list(mtable[mtable['cohort']==negclass]['donor'].unique()))
    donors_l = [donors_pos, donors_neg]
    genes = set(gframe[gframe['geneset']==geneset]['gene'].unique())
    for membership, donors in zip(membership_l, donors_l):
        for donor in donors:
            for gene in sorted(list(genes)):
                cig = f"{membership}:{donor}:{gene}"
                idents.append(cig)
    unstacked = pd.DataFrame(data=0, index=idents, columns=list(cmapper.keys()))

    # populate dataframe
    dfslice = mtable[mtable['gene'].isin(genes)].copy()
    for i, rec in dfslice.iterrows():
        cig = f"{rec.cohort}:{rec.donor}:{rec.gene}"
        unstacked.loc[cig, rec.vclass] = 1
    assert unstacked.isna().sum().sum() == 0

    # resetting df to original format
    unstacked = unstacked.reset_index()
    unstacked['cohort'] = unstacked['index'].apply(lambda x: x.split(':')[0])
    unstacked['donor'] = unstacked['index'].apply(lambda x: x.split(':')[1])
    unstacked['gene'] = unstacked['index'].apply(lambda x: x.split(':')[2])
    unstacked = unstacked.sort_values('cohort')
    unstacked = unstacked.drop('index', axis=1)
    return unstacked

def render_oncoprint(unstacked: pd.DataFrame, col_split: pd.DataFrame, tmb_frame: pd.DataFrame, posclass: str, negclass: str, cmapper: dict, filepath: str) -> None:
    col_ha = pycomp.HeatmapAnnotation(
        axis=1,
        Cohort=pycomp.anno_simple(col_split['cohort'], add_text=True, height=4, cmap='Accent', legend=False),
        Total_Mutation_Counts=pycomp.anno_barplot(
            tmb_frame, 
            colors=[cmapper[x] for x in tmb_frame.columns.to_list()],
            legend=False, height=30, linewidth=0.1
        ),
        verbose=0, 
        label_side='right', 
        hgap=0,
        label_kws={'horizontalalignment': 'left','visible':True}
    )
    
    # main figure
    mtype_order = []
    mtype_colors = []
    for mtype, color in cmapper.items():
        mtype_order.append(mtype)
        mtype_colors.append(color)

    # plot
    N_COLS = col_split.shape[0]
    N_ROWS = unstacked['gene'].nunique()
    WIDTH = min(N_COLS/30 + 2, 40)
    HEIGHT = N_ROWS/4 + 1
    plt.figure(figsize=(WIDTH, HEIGHT), dpi=200)
    op=pycomp.oncoPrintPlotter(
        data=unstacked,
        x='donor',
        y='gene',
        values=mtype_order,
        colors=mtype_colors,
        subplot_gap=4,
        label='Alteration',
        top_annotation=col_ha,
        col_split=col_split['cohort'],
        col_split_order=[posclass, negclass],
        col_split_gap=4,
        row_split_gap=4,
        row_cluster=True,
        col_cluster=True,
        legend_hpad=0,
        show_rownames=True,
        show_colnames=False,
        remove_empty_columns=True
        # remove_empty_columns=False
    )
    plt.savefig(filepath, bbox_inches='tight')
    plt.close()


def load_cmdline_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Summarises contrast results.')
    parser.add_argument('--posclass', type=str, required=True, help='Which group in samplesheet is considered the positive class.')
    parser.add_argument('--genesets', type=str, required=True, help='Path to genesets file (after harmonisation).')
    parser.add_argument('--selection', type=str, required=True, help='Path to genesets selection file (tsv).')
    parser.add_argument('--posmuts', type=str, required=True, help='Path to positive class mutations (tsv).')
    parser.add_argument('--negmuts', type=str, required=True, help='Path to negative class mutations (tsv).')
    parser.add_argument('--diffmut', type=str, required=True, help='Path to differential group test results (tsv).')
    parser.add_argument('--posenrich', type=str, required=True, help='Path to positive class background enrichment results (tsv).')
    parser.add_argument('--negenrich', type=str, required=True, help='Path to negative class background enrichment results (tsv).')
    parser.add_argument('--max-plots', type=int, default=20, help='How many genesets to plot (top <--max-plots>).')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    main()