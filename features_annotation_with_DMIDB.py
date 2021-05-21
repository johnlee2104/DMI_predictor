# This python script contains all the function that calculates and adds the features to PRS and RRS tables.
# Perhaps for sequence i can retrieve them from the previous Protein objects by inheriting from DMIDB.py
import pandas as pd
import time
from DMIDB import *
import DMIDB

""" This script takes a directory and searches for feature files in that directory and read respective features into PRS and RRS.
Original PRS and RRS file should contain only proteins, SLiM and domain annotations (including domain evalues)
and without SLiM features. Using the annotate_slim_features will annotate SLiM features into the PRS and RRS. """
""" dir_name specifies the path where protein features are stored in separate subdirectories."""

def annotate_slim_domain_features_on_dataset(input_list):
    for input_file in input_list:
        df= pd.read_csv(input_file, sep= '\t', index_col= 0)
        new_cols= ['IUPredShort', 'Anchor', 'DomainOverlap', 'qfo_RLC', 'qfo_RLCvar', 'vertebrates_RLC', 'vertebrates_RLCvar', 'mammalia_RLC', 'mammalia_RLCvar', 'metazoa_RLC', 'metazoa_RLCvar', 'DomainEnrichment_pvalue', 'DomainEnrichment_zscore', 'TotalNetworkDegree', 'vertex_with_domain_in_real_network', 'DomainFreqbyProtein1', 'DomainFreqbyProtein2', 'DomainFreqinProteome1', 'DomainFreqinProteome2', 'DMISource']
        for col in new_cols:
            if col not in df.columns:
                df[col]= float('NaN')
        for ind, row in df.iterrows():
            if (pd.notna(row['qfo_RLC'])) & (type(row['qfo_RLC']) == float): # this is when I rerun the code because some RLC values are not calculated due to server not responding
                continue
            if ind % 50 == 0:
                print('Sleeping for 10 seconds to not flood the SLiM server...')
                time.sleep(10)
            slim_protein= row['interactorElm']
            slim_id= row['Accession']
            start, end= row['ElmMatch'].split('-')
            start= int(start)
            end= int(end)
            pattern= row['Pattern'].replace('"', '')
            slim_match_inst= SLiMMatch(InterfaceHandling.dmi_types_dict[slim_id], InterfaceHandling.slim_types_dict[slim_id], InterfaceHandling.proteins_dict[slim_protein], start, end, pattern)
            if pd.notna(row['DomainID2']):
                slim_match_inst.get_slim_match_features(domain_type_list= [row['DomainID1'], row['DomainID2']])
            else:
                slim_match_inst.get_slim_match_features(domain_type_list= [row['DomainID1']])
            features_dict= slim_match_inst.__dict__
            for feature in list(features_dict)[6:]:
                if features_dict[feature] != None:
                    df.loc[ind, feature]= features_dict[feature]
            df.loc[ind, 'TotalNetworkDegree']= InterfaceHandling.proteins_dict[slim_protein].network_degree
            # df.loc[ind, 'DomainFreqbyProtein1']=  InterfaceHandling.domain_types_dict[row['DomainID1']].DomainFreqbyProtein
            # df.loc[ind, 'DomainFreqinProteome1']=  InterfaceHandling.domain_types_dict[row['DomainID1']].DomainFreqinProteome
            # if pd.notna(row['DomainID2']):
            #     df.loc[ind, 'DomainFreqbyProtein2']=  InterfaceHandling.domain_types_dict[row['DomainID2']].DomainFreqbyProtein
            #     df.loc[ind, 'DomainFreqinProteome2']=  InterfaceHandling.domain_types_dict[row['DomainID2']].DomainFreqinProteome

        columns= ['Accession', 'Elm', 'Regex', 'Pattern', 'Probability', 'interactorElm', 'ElmMatch', 'IUPredShort', 'Anchor', 'DomainOverlap', 'qfo_RLC', 'qfo_RLCvar', 'vertebrates_RLC', 'vertebrates_RLCvar', 'mammalia_RLC', 'mammalia_RLCvar', 'metazoa_RLC', 'metazoa_RLCvar', 'DomainEnrichment_pvalue', 'DomainEnrichment_zscore', 'TotalNetworkDegree', 'vertex_with_domain_in_real_network', 'interactorDomain', 'DomainID1', 'DomainMatch1', 'DomainMatchEvalue1', 'DomainFreqbyProtein1', 'DomainFreqinProteome1', 'DomainID2', 'DomainMatch2', 'DomainMatchEvalue2', 'DomainFreqbyProtein2', 'DomainFreqinProteome2', 'DMISource']
        df= df[columns]
        df.to_csv(input_file[:-4] + '_slim_domain_features_annotated.tsv', sep= '\t')
        print(f'New file saved as {input_file[:-4]}_slim_domain_features_annotated.tsv')

if __name__== '__main__':
    prot_path= sys.argv[1]
    PPI_file= sys.argv[2]
    slim_type_file= sys.argv[3]
    dmi_type_file= sys.argv[4]
    smart_domain_types_file= sys.argv[5]
    pfam_domain_types_file= sys.argv[6]
    smart_domain_matches_json_file= sys.argv[7]
    pfam_domain_matches_json_file= sys.argv[8]
    features_path= sys.argv[9]
    network_path= sys.argv[10]
    input_list= list(sys.argv[11:])


    InterfaceHandling= DMIDB.InterfaceHandling(prot_path, slim_type_file, dmi_type_file, smart_domain_types_file, pfam_domain_types_file, smart_domain_matches_json_file, pfam_domain_matches_json_file, features_path, network_path= network_path)

    InterfaceHandling.read_in_proteins()
    InterfaceHandling.read_in_slim_types()
    InterfaceHandling.read_in_DMI_types()
    InterfaceHandling.read_in_domain_types()
    InterfaceHandling.read_in_domain_matches()
    InterfaceHandling.read_in_networks()
    InterfaceHandling.read_in_features_all_proteins()
    annotate_slim_domain_features_on_dataset(input_list)

    # python3 features_annotation_with_DMIDB.py ~/Coding/Python/DMI/protein_sequences_and_features/human_protein_sequences PRS_IntAct_union_known_PPIs.txt ../elm_classes_20210222.tsv ../elm_interaction_domains_complete_20210222.tsv ../domain_stuffs/all_smart_domains_with_frequency.txt ../domain_stuffs/all_pfam_domains_with_frequency.txt ../domain_stuffs/interpro_9606_smart_matches_20210122.json ../domain_stuffs/interpro_9606_pfam_matches_20210122.json ../protein_sequences_and_features/human_protein_sequences_features ../protein_sequences_and_features/human_protein_sequences_features/Protein_networks_PRS_filtered PRS_v3_only_human_with_pattern_alt_iso_swapped_removed_20210413.tsv

    # python3 features_annotation_with_DMIDB.py ~/Coding/Python/DMI/protein_sequences_and_features/human_protein_sequences PRS_IntAct_union_known_PPIs.txt ../elm_classes_20210222.tsv ../elm_interaction_domains_complete_20210222.tsv ../domain_stuffs/all_smart_domains_with_frequency.txt ../domain_stuffs/all_pfam_domains_with_frequency.txt ../domain_stuffs/interpro_9606_smart_matches_20210122.json ../domain_stuffs/interpro_9606_pfam_matches_20210122.json ../protein_sequences_and_features/human_protein_sequences_features ../protein_sequences_and_features/human_protein_sequences_features/Protein_networks_PRS_filtered RRSv1_1_20210427.tsv RRSv1_2_20210427.tsv RRSv1_3_20210427.tsv RRSv2_1_20210428.tsv RRSv2_2_20210428.tsv RRSv2_3_20210428.tsv RRSv3_1_20210428.tsv RRSv3_2_20210428.tsv RRSv3_3_20210428.tsv RRSv4_1_20210428.tsv RRSv4_2_20210428.tsv RRSv4_3_20210428.tsv
