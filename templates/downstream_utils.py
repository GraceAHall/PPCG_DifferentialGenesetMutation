
import numpy as np 
import pandas as pd
 

#####################
### hypermutators ###
#####################

def filter_hypermutators(df: pd.DataFrame, zscore_thresh: float):
    counts = df.drop_duplicates(subset=['donor', 'vclass', 'gene'])['donor'].value_counts()
  
    # Modified Z-Score (MAD)
    median = counts.median()
    mad = np.median(np.abs(counts - median))
    
    # Modified Z-score
    modified_z_scores = 0.6745 * (counts - median) / mad
    
    # Summary
    results = counts.to_frame()
    results['modified_z_scores'] = modified_z_scores
    results['outlier'] = results['modified_z_scores']>zscore_thresh
    hypermutators = set(results[results['outlier']==True].index.to_list())
    df_filt = df[~df['donor'].isin(hypermutators)].copy()

    return df_filt, results