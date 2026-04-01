
import pandas as pd 
import seaborn as sns 
import matplotlib.pyplot as plt 


def summarise_basic_info(df: pd.DataFrame) -> None:
    print('\n--- Basic Info ---\n')
    print()
    sframe = pd.DataFrame(index=df['cohort'].unique())
    sframe['Donors'] = df.groupby('cohort')['donor'].nunique()
    sframe['Genes'] = df.groupby('cohort')['gene'].nunique()
    sframe['Variants'] = df.groupby('cohort')['ID'].nunique()
    print(sframe)

def summarise_donors(df: pd.DataFrame) -> None:
    print('\n--- Donors ---')
    for clsmem in df['cohort'].unique():
        print()
        print(f'[{clsmem}]')
        # dfslice = df[df['cohort']==clsmem].drop_duplicates('ID')
        dfslice = df[df['cohort']==clsmem].drop_duplicates(['donor', 'vclass', 'gene'])
        counts = dfslice.groupby('donor')['vclass'].value_counts().unstack().fillna(0).astype(int)
        counts['total'] = dfslice.groupby('donor')['gene'].nunique()
        counts = counts.sort_values('total', ascending=False)
        print(counts.head(10))

def summarise_vclasses(df: pd.DataFrame) -> None:
    print('\n--- Variant Class ---')
    vclasses = sorted(list(df['vclass'].unique()))

    for clsmem in df['cohort'].unique():
        print()
        print(f'[{clsmem}]')
        dfslice = df[df['cohort']==clsmem].drop_duplicates('ID')
        counts = pd.DataFrame(index=vclasses)
        counts['vars'] = dfslice['vclass'].value_counts()
        counts['genes'] = dfslice.groupby('vclass')['gene'].nunique()
        counts['donors'] = dfslice.groupby('vclass')['donor'].nunique()
        
        for vclass in vclasses:
            counts.loc[vclass, 'med. vars p/donor'] = round(dfslice[dfslice['vclass']==vclass]['donor'].value_counts().median(), 0)
        counts = counts.astype(int)
        print(counts)

def summarise_annotations(df: pd.DataFrame) -> None:
    print('\n--- Variant Annotations ---')
    annotations = sorted(list(df['annotation'].unique()))

    for clsmem in df['cohort'].unique():
        print()
        print(f'[{clsmem}]')
        dfslice = df[df['cohort']==clsmem].drop_duplicates('ID')
        counts = pd.DataFrame(index=annotations)
        counts['vars'] = dfslice['annotation'].value_counts()
        counts['genes'] = dfslice.groupby('annotation')['gene'].nunique()
        counts['donors'] = dfslice.groupby('annotation')['donor'].nunique()
        
        for annot in annotations:
            counts.loc[annot, 'med. vars p/donor'] = round(dfslice[dfslice['annotation']==annot]['donor'].value_counts().median(), 0)
        counts = counts.astype(int)
        print(counts)

def top_genes(df: pd.DataFrame) -> None:
    print('\n--- Gene Summary (stratified by variant class) ---')
    for clsmem in df['cohort'].unique():
        print()
        print(f'[{clsmem}]')
        dfslice = df[df['cohort']==clsmem]
        ndonors = dfslice['donor'].nunique()
        all_genes = sorted(list(dfslice['gene'].unique()))
        counts = pd.DataFrame(index=all_genes)
        counts['sv'] = dfslice[dfslice['vclass']=='SV'].groupby('gene')['donor'].nunique()
        counts['snvs'] = dfslice[dfslice['vclass']=='SNV'].groupby('gene')['donor'].nunique()
        counts['indels'] = dfslice[dfslice['vclass']=='INDEL'].groupby('gene')['donor'].nunique()
        counts['cna'] = dfslice[dfslice['vclass']=='CNA'].groupby('gene')['donor'].nunique()
        counts['total'] = dfslice.groupby('gene')['donor'].nunique()
        counts = counts.fillna(0).astype(int)
        counts['donors prop.'] = counts['total'] / ndonors
        counts['donors prop.'] = counts['donors prop.'].apply(lambda x: f"{x*100:.1f} %")
        counts = counts.sort_values(by='total', ascending=False)
        print(counts.head(20))

def plot_distribution(df: pd.DataFrame, filepath: str) -> None:
    title = 'Number of Genes Mutated per Donor'
    counts = df.groupby('donor')['gene'].nunique().to_frame()
    counts['cohort'] = df.drop_duplicates('donor').set_index('donor')['cohort'].to_dict()
    counts = counts.reset_index()
    sns.histplot(counts, x='gene', hue='cohort', bins=10, multiple='dodge', stat='proportion', common_norm=False)
    plt.title(title)
    plt.savefig(filepath)
    plt.close()
