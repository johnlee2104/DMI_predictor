# This python script contains all the function that calculates and adds the features to PRS and RRS tables.
# Perhaps for sequence i can retrieve them from the previous Protein objects by inheriting from DMIDB.py
import pandas as pd
import sys
import glob
""" This script takes a directory and searches for feature files in that directory and read respective features into PRS and RRS.
Original PRS and RRS file should contain only proteins, SLiM and domain annotations (including domain evalues)
and without SLiM features. Using the annotate_slim_features will annotate SLiM features into the PRS and RRS. """

def annotate_slim_features_on_dataset(dir_name, input_file_tsv):
    df= pd.read_csv(input_file_tsv, sep= '\t', index_col= 0)
    df['IUPredLong']= float('NaN')
    df['IUPredShort']= float('NaN')
    df['Anchor']= float('NaN')
    for ind, row in df.iterrows():
        protein_id= row['interactorElm']
        start, end= row['ElmMatch'].split('-')
        start= int(start)
        end= int(end)
        IUPredLong_score= []
        IUPredShort_score= []
        Anchor_score= []
        for file_name in glob.glob(dir_name + '/**/' + protein_id + '*.txt', recursive= True):
            with open(file_name, 'r') as f:
                lines= [line.strip() for line in f.readlines()]
                for line in lines[1:]:
                    if 'long' in file_name:
                        IUPredLong_score.append(float(line.split('\t')[2]))
                    elif 'short' in file_name:
                        IUPredShort_score.append(float(line.split('\t')[2]))
                    elif 'anchor' in file_name:
                        Anchor_score.append(float(line.split('\t')[2]))
        df.loc[ind, ['IUPredLong']]= format(sum(IUPredLong_score[start-1:end])/(end-start+1), '.4f')
        df.loc[ind, ['IUPredShort']]= format(sum(IUPredShort_score[start-1:end])/(end-start+1), '.4f')
        df.loc[ind, ['Anchor']]= format(sum(Anchor_score[start-1:end])/(end-start+1), '.4f')
        print(f'{protein_id}\'s SLiM features annotated.')
    df.to_csv(input_file_tsv[:-4] + '_slim_features_annotated.tsv', sep= '\t')

if __name__== '__main__':
    dir_name= sys.argv[1]
    input_tsv_file= sys.argv[2]
    annotate_slim_features_on_dataset(dir_name, input_tsv_file)
