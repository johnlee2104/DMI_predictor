# This script trains the final random forest model for deployment using ALL the data in the PRS and the selected RRS version.
# Author: Chop Yan Lee

import sys
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.tree import plot_tree
from joblib import dump

all_features= ['Probability', 'IUPredShort', 'Anchor', 'DomainOverlap', 'qfo_RLC', 'qfo_RLCvar', 'vertebrates_RLC', 'vertebrates_RLCvar', 'mammalia_RLC', 'mammalia_RLCvar', 'metazoa_RLC', 'metazoa_RLCvar', 'DomainEnrichment_pvalue', 'DomainEnrichment_zscore', 'DomainFreqbyProtein1', 'DomainFreqinProteome1']
all_features_renamed= ['Probability', 'IUPredShort', 'Anchor', 'DomainOverlap', 'qfo_RLC', 'qfo_RLCvar', 'vertebrates_RLC', 'vertebrates_RLCvar', 'mammalia_RLC', 'mammalia_RLCvar', 'metazoa_RLC', 'metazoa_RLCvar', 'DomainEnrichment_pvalue', 'DomainEnrichment_zscore', 'DomainFreqbyProtein', 'DomainFreqinProteome']

def fit_final_model_imputer(PRS_path, RRS_path, save_path):
    """
    Fits a random forest classifier and a median imputer to the provided PRS and RRS datasets using ALL their data points (Data points with any NaN value is removed from the dataset). As some DMI requires two domains from the domain protein to bind to a SLiM, the domain features of the two domains are averaged to produce one domain feature for training. The concatenated PRS and RRS with NaN rows removed is saved as a tsv file.

    Args:
        PRS_path (str): Absolute path to the PRS dataset
        RRS_path (str): Absolute path to the RRS dataset
        save_path (str): Absolute path to where the concatenated dataframe should be saved

    Returns:
        rf (sklearn.ensemble.RandomForestClassifier): The fitted random forest model
        median_imp (sklearn.impute.SimpleImputer): The fitted median imputer
    """
    PRS= pd.read_csv(PRS_path, sep= '\t', index_col= 0)
    PRS['label']= 1
    RRS= pd.read_csv(RRS_path, sep= '\t', index_col= 0)
    RRS['label']= 0
    for df in [PRS, RRS]:
        for ind, row in df.iterrows():
            if pd.notna(row['DomainFreqbyProtein2']):
                df.loc[ind, 'DomainFreqbyProtein1'] = np.mean([row['DomainFreqbyProtein1'], row['DomainFreqbyProtein2']])
                df.loc[ind, 'DomainFreqinProteome1'] = np.mean([row['DomainFreqinProteome1'], row['DomainFreqinProteome2']])
    df= pd.concat([PRS, RRS], ignore_index= True)
    df.dropna(subset= all_features, inplace= True)
    df.rename(columns= {'DomainFreqbyProtein1': 'DomainFreqbyProtein', 'DomainFreqinProteome1': 'DomainFreqinProteome'}, inplace= True)
    X= df[all_features_renamed].copy()
    y= df['label']

    rf= RandomForestClassifier(n_estimators= 1000, oob_score= True, verbose= True, n_jobs= 1) # Tested the model again on 11.03.2024 and found that it only works if n_jobs is set to 1. -1 makes use of all CPU cores but somehow it does not work on my computer.
    rf.fit(X, y)

    median_imp= SimpleImputer(strategy= 'median')
    median_imp.fit(X)

    df.to_csv(save_path + f'final_PRS_{RRS_version}_to_fit_model.tsv', sep= '\t')

    return rf, median_imp

def save_model_imputer(model, imputer, save_path):
    """
    Save the fitted model and imputer as joblib files.

    Args:
        model (sklearn.ensemble.RandomForestClassifier): The fitted random forest model
        imputer (sklearn.impute.SimpleImputer): The fitted median imputer
        save_path (str): Absolute path where the model and imputer should be saved
    """
    with open(save_path + f'final_RF_model_with_{RRS_version}.joblib', 'wb') as f:
        dump(model, f)

    with open(save_path + f'final_median_imputer_with_{RRS_version}.joblib', 'wb') as f:
        dump(imputer, f)

if __name__ == '__main__':
    PRS= sys.argv[1]
    RRS= sys.argv[2]
    RRS_version= RRS.split('/')[-1][:7]
    save_path= sys.argv[3]

    rf, median_imp= fit_final_model_imputer(PRS, RRS, save_path)
    save_model_imputer(rf, median_imp, save_path)

# python3 fitting_final_RF_imputer.py ../PRS/PRS_v3_only_human_with_pattern_alt_iso_swapped_removed_20210413_slim_domain_features_annotated.tsv ../RRS/RRSv4/RRSv4_3_20210428_slim_domain_features_annotated.tsv ~/Coding/Python/DMI/
