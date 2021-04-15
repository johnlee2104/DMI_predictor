# This python script contains all the function that calculates and adds the features to PRS and RRS tables.
# Perhaps for sequence i can retrieve them from the previous Protein objects by inheriting from DMIDB.py
import pandas as pd
import sys, glob, os, requests, json
import scipy.special as sc
import numpy as np
from DMIDB import *

""" This script takes a directory and searches for feature files in that directory and read respective features into PRS and RRS.
Original PRS and RRS file should contain only proteins, SLiM and domain annotations (including domain evalues)
and without SLiM features. Using the annotate_slim_features will annotate SLiM features into the PRS and RRS. """
""" dir_name specifies the path where protein features are stored in separate subdirectories."""

# global variable
glob_input_file= input_file_tsv

def annotate_IUP_domainoverlap_domainfreqs_on_dataset(slim_features_dir_name, domain_features_dir_name, input_file_tsv):
    domain_freq_dict= {}
    for file_name in glob.glob(domain_features_dir_name + '/*_with_frequency.txt'):
        print(file_name)
        with open(file_name, 'r') as f:
            lines= [line.strip() for line in f.readlines()]
        for line in lines[2:]:
            tabs= line.split('\t')
            if len(tabs)== 5:
                domain_freq_dict[tabs[1]]= []
                domain_freq_dict[tabs[1]].append(tabs[3])
                domain_freq_dict[tabs[1]].append(tabs[4])
    df= pd.read_csv(input_file_tsv, sep= '\t', index_col= 0)
    new_cols= ['IUPredLong', 'IUPredShort', 'Anchor', 'DomainOverlap', 'DomainFreqbyProtein1', 'DomainFreqbyProtein2', 'DomainFreqinProteome1', 'DomainFreqinProteome2', 'DMISource']
    for col in new_cols:
        if col not in df.columns:
            df[col]= float('NaN')
    for ind, row in df.iterrows():
        protein_id= row['interactorElm']
        regex= row['Regex']
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
    columns= ['Accession', 'Elm', 'Regex', 'Pattern', 'Probability', 'interactorElm', 'ElmMatch', 'IUPredLong', 'IUPredShort', 'Anchor', 'DomainOverlap', 'interactorDomain', 'DomainID1', 'DomainMatch1', 'DomainMatchEvalue1', 'DomainFreqbyProtein1', 'DomainFreqinProteome1', 'DomainID2', 'DomainMatch2', 'DomainMatchEvalue2', 'DomainFreqbyProtein2', 'DomainFreqinProteome2', 'DMISource']
    df= df[columns]
    target= input_file_tsv[:-4] + '_IUP_domainoverlap_domainfreqs.tsv'
    df.to_csv(target, sep= '\t')
    print(f'New file saved as {target}')
    # assign the new file name to global variable so the new file will be used for cons score annotation
    global glob_input_file
    glob_input_file= target

def annotate_slim_cons_scores_on_dataset(cons_score_path, input_file_tsv= None):
    if input_file_tsv == None:
        input_file_tsv= glob_input_file
    else:
        input_file_tsv= input_file_tsv
    base_url= 'http://slim.icr.ac.uk/restapi/functions/defined_positions?'
    df= pd.read_csv(input_file_tsv, sep= '\t', index_col= 0)
    new_cols= ['qfo_RLC', 'qfo_RLCvar', 'vertebrates_RLC', 'vertebrates_RLCvar', 'mammalia_RLC', 'mammalia_RLCvar', 'metazoa_RLC', 'metazoa_RLCvar', 'DMISource']
    for col in new_cols:
        if col not in df.columns:
            df[col]= float('NaN')
    for ind, row in df.iterrows():
        if row['qfo_RLC'] != 'Server still not responding after 61 secs timeout':
            continue
        protein_id= row['interactorElm']
        regex= row['Regex']
        pattern= row['Pattern']
        payload= {'motif':regex, 'sequence':pattern}
        start, end= row['ElmMatch'].split('-')
        start= int(start)
        end= int(end)
        if protein_id + '_con.json' in os.listdir(cons_score_path + '/'):
            try:
                r= requests.get(base_url, params= payload, timeout= 61)
                if r.status_code == requests.codes.ok:
                    response= r.json()
                    defined_positions= [start + (ind - 1) for ind in response['indexes']]
                    # ind - 1 because the pattern sent is always with flanking residues of the regex match
                with open(cons_score_path + '/' + protein_id + '_con.json', 'r') as f:
                    data= json.load(f) # pos annotated starting from 1
                    for result in data['Conservation']:
                        for cons_type in result:
                            defined_positions_cons_scores= []
                            for pos, score in result[cons_type].items():
                                if int(pos) in defined_positions:
                                    defined_positions_cons_scores.append(score)
                            if any(defined_positions_cons_scores):
                                pmotif= np.product(defined_positions_cons_scores)
                                lnpmotif= -np.log(pmotif)
                                sigmotif= sc.gammaincc(len(defined_positions_cons_scores), lnpmotif)
                                meanRLCprob= np.mean(defined_positions_cons_scores)
                                varRLCprob= sum([abs(x-meanRLCprob) for x in defined_positions_cons_scores])/len(defined_positions_cons_scores)
                                df.loc[ind, [cons_type + '_RLC']]= sigmotif
                                df.loc[ind, [cons_type + '_RLCvar']]= varRLCprob
                                print(f'{sigmotif} & {varRLCprob} annotated for {cons_type} of {protein_id}.')
                            else:
                                df.loc[ind, [cons_type + '_RLC']]= 'Problem with regex/motif'
                                df.loc[ind, [cons_type + '_RLCvar']]= 'Problem with regex/motif'
            except:
                df.loc[ind, 'qfo_RLC']= 'Server still not responding after 61 secs timeout'
        else:
            df.loc[ind, 'qfo_RLC']= 'No conservation score found'
    columns= ['Accession', 'Elm', 'Regex', 'Pattern', 'Probability', 'interactorElm', 'ElmMatch', 'IUPredLong', 'IUPredShort', 'Anchor', 'DomainOverlap', 'qfo_RLC', 'qfo_RLCvar', 'vertebrates_RLC', 'vertebrates_RLCvar', 'mammalia_RLC', 'mammalia_RLCvar', 'metazoa_RLC', 'metazoa_RLCvar', 'interactorDomain', 'DomainID1', 'DomainMatch1', 'DomainMatchEvalue1', 'DomainFreqbyProtein1', 'DomainFreqinProteome1', 'DomainID2', 'DomainMatch2', 'DomainMatchEvalue2', 'DomainFreqbyProtein2', 'DomainFreqinProteome2', 'DMISource']
    df= df[columns]
    target= input_file_tsv[:-4] + '_cons_scores.tsv'
    df.to_csv(target, sep= '\t')
    print(f'New file saved as {target}')
    # assign the new file name to global variable so the new file will be used for cons score annotation
    global glob_input_file
    glob_input_file= target

def annotate_network_features_on_dataset(network_path, input_file_tsv= None):
    if input_file_tsv == None:
        input_file_tsv= glob_input_file
    else:
        input_file_tsv= input_file_tsv


if __name__== '__main__':
    prot_path= sys.argv[1]
    slim_type_file= sys.argv[2]
    dmi_type_file= sys.argv[3]
    smart_domain_types_file= sys.argv[4]
    pfam_domain_types_file= sys.argv[5]
    smart_domain_matches_json_file= sys.argv[6]
    pfam_domain_matches_json_file= sys.argv[7]
    features_path= sys.argv[8]
    slim_features_dir_name= sys.argv[1]
    domain_features_dir_name= sys.argv[2]
    cons_score_path= sys.argv[3]
    input_file_tsv= sys.argv[4]

    InterfaceHandling= DMIDB.InterfaceHandling(slim_type_file, dmi_type_file, smart_domain_types_file, pfam_domain_types_file, smart_domain_matches_json_file, pfam_domain_matches_json_file, features_path)
    InterfaceHandling.read_in_slim_types()
    InterfaceHandling.read_in_DMI_types()
    InterfaceHandling.read_in_domain_types()
    annotate_slim_domain_features_on_dataset(slim_features_dir_name, domain_features_dir_name, cons_score_path, input_file_tsv)

    # python3 features_annotation.py ../protein_sequences_and_features/PRS_v3_RRS_v1_sequences_features ../domain_stuffs ../protein_sequences_and_features/human_protein_sequences_features/conservation_scores RRSv1_20210315.tsv
