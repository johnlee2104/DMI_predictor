# This python script contains all the function that calculates and adds the features to PRS and RRS tables.
# Perhaps for sequence i can retrieve them from the previous Protein objects by inheriting from DMIDB.py
import pandas as pd
import sys, glob

""" This script takes a directory and searches for feature files in that directory and read respective features into PRS and RRS.
Original PRS and RRS file should contain only proteins, SLiM and domain annotations (including domain evalues)
and without SLiM features. Using the annotate_slim_features will annotate SLiM features into the PRS and RRS. """
""" dir_name specifies the path where protein features are stored in separate subdirectories."""
def annotate_slim_domain_features_on_dataset(slim_features_dir_name, domain_features_dir_name, input_file_tsv):
    domain_freq_dict= {}
    for file_name in glob.glob(domain_features_dir_name + '/*_with_frequency.txt'):
        with open(file_name, 'r') as f:
            lines= [line.strip() for line in f.readlines()]
        for line in lines[2:]:
            tabs= line.split('\t')
            if len(tabs)== 5:
                domain_freq_dict[tabs[1]]= []
                domain_freq_dict[tabs[1]].append(tabs[3])
                domain_freq_dict[tabs[1]].append(tabs[4])
    print(domain_freq_dict['PF03568'])
    df= pd.read_csv(input_file_tsv, sep= '\t', index_col= 0)
    df['IUPredLong']= float('NaN')
    df['IUPredShort']= float('NaN')
    df['Anchor']= float('NaN')
    df['DomainOverlap']= float('NaN')
    df['DomainFreqbyProtein1']= float('NaN')
    df['DomainFreqinProteome1']= float('NaN')
    df['DomainFreqbyProtein2']= float('NaN')
    df['DomainFreqinProteome2']= float('NaN')
    if 'DMISource' not in df.columns:
        df['DMISource']= float('NaN')
    for ind, row in df.iterrows():
        protein_id= row['interactorElm']
        start, end= row['ElmMatch'].split('-')
        domain_id1= row['DomainID1']
        domain_id2= ''
        if pd.isna(row['DomainID2'])== False:
            domain_id2= row['DomainID2']
        start= int(start)
        end= int(end)
        IUPredLong_score= []
        IUPredShort_score= []
        Anchor_score= []
        DomainOverlap_score= []
        with open(slim_features_dir_name + '/IUPred_long/' + protein_id + '_iupredlong.txt', 'r') as f:
            lines= [line.strip() for line in f.readlines()]
            for line in lines[1:]:
                IUPredLong_score.append(float(line.split('\t')[2]))
        with open(slim_features_dir_name + '/IUPred_short/' + protein_id + '_iupredshort.txt', 'r') as f:
            lines= [line.strip() for line in f.readlines()]
            for line in lines[1:]:
                IUPredShort_score.append(float(line.split('\t')[2]))
        with open(slim_features_dir_name + '/Anchor/' + protein_id + '_anchor.txt', 'r') as f:
            lines= [line.strip() for line in f.readlines()]
            for line in lines[1:]:
                Anchor_score.append(float(line.split('\t')[2]))
        with open(slim_features_dir_name + '/Domain_overlap/' + protein_id + '_domain_overlap.txt', 'r') as f:
            lines= [line.strip() for line in f.readlines()]
            for line in lines[1:]:
                DomainOverlap_score.append(float(line.split('\t')[2]))
        df.loc[ind, ['IUPredLong']]= format(sum(IUPredLong_score[start-1:end])/(end-start+1), '.4f')
        df.loc[ind, ['IUPredShort']]= format(sum(IUPredShort_score[start-1:end])/(end-start+1), '.4f')
        df.loc[ind, ['Anchor']]= format(sum(Anchor_score[start-1:end])/(end-start+1), '.4f')
        df.loc[ind, ['DomainOverlap']]= format(sum(DomainOverlap_score[start-1:end])/(end-start+1), '.4f')
        print(f'{protein_id}\'s SLiM features annotated.')
        df.loc[ind, ['DomainFreqbyProtein1']]= domain_freq_dict[domain_id1][0]
        df.loc[ind, ['DomainFreqinProteome1']]= domain_freq_dict[domain_id1][1]
        print(f'{domain_id1}\'s domain features annotated.')
        if pd.isna(row['DomainID2'])== False:
            df.loc[ind, ['DomainFreqbyProtein2']]= domain_freq_dict[domain_id2][0]
            df.loc[ind, ['DomainFreqinProteome2']]= domain_freq_dict[domain_id2][1]
            print(f'{domain_id2}\'s domain features annotated.')
    columns= ['Accession', 'Elm', 'Regex', 'Probability', 'interactorElm', 'ElmMatch', 'IUPredLong', 'IUPredShort', 'Anchor', 'DomainOverlap',
    'interactorDomain', 'DomainID1', 'DomainMatch1', 'DomainMatchEvalue1', 'DomainFreqbyProtein1', 'DomainFreqinProteome1', 'DomainID2',
    'DomainMatch2', 'DomainMatchEvalue2', 'DomainFreqbyProtein2', 'DomainFreqinProteome2', 'DMISource']
    df= df[columns]
    df.to_csv(input_file_tsv[:-4] + '_slim_domain_features_annotated.tsv', sep= '\t')
    print(f'New file saved as {input_file_tsv[:-4]}_slim_domain_features_annotated.tsv')

if __name__== '__main__':
    slim_features_dir_name= sys.argv[1]
    domain_features_dir_name= sys.argv[2]
    input_file_tsv= sys.argv[3]
    annotate_slim_domain_features_on_dataset(slim_features_dir_name, domain_features_dir_name, input_file_tsv)
