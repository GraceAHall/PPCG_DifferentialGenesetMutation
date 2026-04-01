import pandas as pd 


def load_hframe(filepath: str) -> pd.DataFrame:
    sheet = pd.read_csv(filepath, sep='\t', header=0)
    sheet['tissue'] = sheet['tissue'].replace('Recurrence', 'Metastasis')
    sheet = sheet.sort_values('sample')

    # annotating which files are present
    vclasses = ['SNV', 'INDEL', 'CNA', 'SV']
    sheet['Files'] = sheet['Missing_Files'].apply(format_files)
    for vclass in vclasses:
        sheet[vclass] = sheet['Files'].apply(lambda x: vclass in x)
    sheet = sheet.drop(['Files'], axis=1)

    data = []
    cohort_LUT = sheet.drop_duplicates('donor').set_index('donor')['cohort'].to_dict()
    for vclass in vclasses:
        for donor in sheet['donor'].unique():
            label = f"{donor}|{vclass}"
            cohort = cohort_LUT[donor]
            temp = sheet[(sheet['donor']==donor) & (sheet[vclass]==True)]
            primary = temp[temp['tissue']=='Primary']['sample'].nunique()
            secondary = temp[temp['tissue']=='Metastasis']['sample'].nunique()
            data.append((label, donor, vclass, cohort, primary, secondary))
    hframe = pd.DataFrame(data, columns=['label', 'donor', 'vclass', 'cohort', 'primary', 'secondary'])
    hframe['handling'] = hframe.apply(get_handling, axis=1)
    return hframe

def format_files(missing_text: str) -> set[str]:
    expected = set(['SNV', 'INDEL', 'CNA', 'SV'])
    if missing_text == '.':
        return expected
    missing = set(missing_text.split(','))
    return expected-missing

def get_handling(row: pd.Series) -> str:
    if row['primary'] == 0 and row['secondary'] == 0:
        return 'ignore'
    elif row['primary'] == 1 and row['secondary'] == 0:
        return 'single_primary'
    elif row['primary'] >= 2 and row['secondary'] == 0:
        return 'multi_primary'
    elif row['primary'] == 0 and row['secondary'] == 1:
        return 'single_met'
    elif row['primary'] == 0 and row['secondary'] >= 2:
        return 'multi_met'
    elif row['primary'] == 1 and row['secondary'] == 1:
        return 'single_primary_single_met'
    elif row['primary'] == 1 and row['secondary'] >= 2:
        return 'single_primary_multi_met'
    raise NotImplementedError(row)
