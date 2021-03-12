""" Takes a list of interaction pairs and find DMI matches between two interacting proteins. """

import sys, json
sys.path.append('/Users/chopyanlee/Coding/Python/DMI/protein_interaction_prediction/')
from DMIDB import *
import DMIDB
import pandas as pd

def find_DMI_matches_in_interaction_file(prot_path, interaction_file, DMIMatch_outfile):
    df= pd.read_csv(interaction_file, sep= '\t')
    file= open(DMIMatch_outfile + '.tsv', 'w')
    file.write('\t'.join(('intx_ID', 'Accession', 'Elm', 'Regex', 'Pattern', 'Probability', 'interactorElm', 'ElmMatch', 'interactorDomain', 'DomainID1', 'DomainMatch1', 'DomainID2', 'DomainMatch2')))
    file.write('\n')
    df['PotentialDMIMatch']= float('NaN')
    df['domains_1_overlap']= float('NaN')
    df['domains_2_overlap']= float('NaN')
    for i, r in df.iterrows():
        protein_1= r['protein_1']
        protein_2= r['protein_2']
        intx_ID= r['intx_ID']
        seq_start_1= int(r['seq_start_1'])
        seq_end_1= int(r['seq_end_1'])
        seq_start_2= int(r['seq_start_2'])
        seq_end_2= int(r['seq_end_2'])
        domains_1= set()
        if pd.isna(r['domains_1']) == False:
            if ';' in r['domains_1']:
                for domain in r['domains_1'].split(';'):
                    domains_1.add(domain)
            else:
                domains_1.add(r['domains_1'])
        domains_2= set()
        if pd.isna(r['domains_2']) == False:
            if ';' in r['domains_2']:
                for domain in r['domains_2'].split(';'):
                    domains_2.add(domain)
            else:
                domains_1.add(r['domains_2'])
        filtered_slim_matches_1= []
        filtered_slim_matches_2= []
        filtered_domain_matches_1= []
        filtered_domain_matches_2= []
        domains_1_overlap= set()
        domains_2_overlap= set()
        InterfaceHandling.proteins_dict= {}
        InterfaceHandling.protein_pairs_dict= {}
        protein_pair= tuple([protein_1, protein_2])
        InterfaceHandling.protein_pairs_dict[protein_pair]= DMIDB.ProteinPair(protein_pair[0], protein_pair[1])
        for protein_id in [protein_1, protein_2]:
            with open(prot_path + '/' + protein_id + '.txt', 'r') as f:
                lines= [line.strip() for line in f.readlines()]
            prot_inst= Protein(protein_id)
            InterfaceHandling.proteins_dict[protein_id]= prot_inst
            InterfaceHandling.proteins_dict[protein_id].sequence= lines[1]
        InterfaceHandling.create_slim_matches_all_proteins()
        for slim_id, slim_matches in InterfaceHandling.proteins_dict[protein_1].slim_matches_dict.items():
            filtered_slim_matches_1= []
            print(f'{protein_1} has {len(InterfaceHandling.proteins_dict[protein_1].slim_matches_dict[slim_id])} {slim_id} match before')
            for slim_match in slim_matches:
                if set(range(slim_match.start, slim_match.end + 1)).issubset(range(seq_start_1,seq_end_1 + 1)):
                    filtered_slim_matches_1.append(slim_match)
            InterfaceHandling.proteins_dict[protein_1].slim_matches_dict[slim_id]= filtered_slim_matches_1
            print(f'{protein_1} has {len(InterfaceHandling.proteins_dict[protein_1].slim_matches_dict[slim_id])} {slim_id} match after')
        for slim_id, slim_matches in InterfaceHandling.proteins_dict[protein_2].slim_matches_dict.items():
            filtered_slim_matches_2= []
            print(f'{protein_2} has {len(InterfaceHandling.proteins_dict[protein_2].slim_matches_dict[slim_id])} {slim_id} match before')
            for slim_match in slim_matches:
                if set(range(slim_match.start, slim_match.end + 1)).issubset(range(seq_start_2,seq_end_2 + 1)):
                    filtered_slim_matches_2.append(slim_match)
            InterfaceHandling.proteins_dict[protein_2].slim_matches_dict[slim_id]= filtered_slim_matches_2
            print(f'{protein_2} has {len(InterfaceHandling.proteins_dict[protein_2].slim_matches_dict[slim_id])} {slim_id} match after')
        for slim_id in list(InterfaceHandling.proteins_dict[protein_1].slim_matches_dict):
            if len(InterfaceHandling.proteins_dict[protein_1].slim_matches_dict[slim_id])== 0:
                del InterfaceHandling.proteins_dict[protein_1].slim_matches_dict[slim_id]
        for slim_id in list(InterfaceHandling.proteins_dict[protein_2].slim_matches_dict):
            if len(InterfaceHandling.proteins_dict[protein_2].slim_matches_dict[slim_id])== 0:
                del InterfaceHandling.proteins_dict[protein_2].slim_matches_dict[slim_id]
        InterfaceHandling.read_in_domain_matches()
        for domain_id, domain_matches in InterfaceHandling.proteins_dict[protein_1].domain_matches_dict.items():
            filtered_domain_matches_1= []
            print(f'{protein_1} has {len(InterfaceHandling.proteins_dict[protein_1].domain_matches_dict[domain_id])} {domain_id} match before')
            for domain_match in domain_matches:
                if (len(set(range(domain_match.start, domain_match.end)).intersection(set(range(seq_start_1, seq_end_1))))/len(range(seq_start_1, seq_end_1)) > 0.8) | (set(range(domain_match.start, domain_match.end)).issubset(range(seq_start_1, seq_end_1 + 1))):
                    filtered_domain_matches_1.append(domain_match)
                    print(domain_match.start, domain_match.end)
                    print(seq_start_1, seq_end_1)
            InterfaceHandling.proteins_dict[protein_1].domain_matches_dict[domain_id]= filtered_domain_matches_1
            print(f'{protein_1} has {len(InterfaceHandling.proteins_dict[protein_1].domain_matches_dict[domain_id])} {domain_id} match after')
        for domain_id, domain_matches in InterfaceHandling.proteins_dict[protein_2].domain_matches_dict.items():
            filtered_domain_matches_2= []
            print(f'{protein_2} has {len(InterfaceHandling.proteins_dict[protein_2].domain_matches_dict[domain_id])} {domain_id} match before')
            for domain_match in domain_matches:
                if (len(set(range(domain_match.start, domain_match.end)).intersection(set(range(seq_start_2, seq_end_2))))/len(range(seq_start_2, seq_end_2)) > 0.8) | (set(range(domain_match.start, domain_match.end)).issubset(range(seq_start_2, seq_end_2 + 1))):
                    filtered_domain_matches_2.append(domain_match)
                    print(domain_match.start, domain_match.end)
                    print(seq_start_2, seq_end_2)
            InterfaceHandling.proteins_dict[protein_2].domain_matches_dict[domain_id]= filtered_domain_matches_2
            print(f'{protein_2} has {len(InterfaceHandling.proteins_dict[protein_2].domain_matches_dict[domain_id])} {domain_id} match after')
        for domain_id in list(InterfaceHandling.proteins_dict[protein_1].domain_matches_dict):
            if len(InterfaceHandling.proteins_dict[protein_1].domain_matches_dict[domain_id])== 0:
                del InterfaceHandling.proteins_dict[protein_1].domain_matches_dict[domain_id]
        for domain_id in list(InterfaceHandling.proteins_dict[protein_2].domain_matches_dict):
            if len(InterfaceHandling.proteins_dict[protein_2].domain_matches_dict[domain_id])== 0:
                del InterfaceHandling.proteins_dict[protein_2].domain_matches_dict[domain_id]
        InterfaceHandling.find_DMI_matches()
        if any(InterfaceHandling.protein_pairs_dict[protein_pair].dmi_matches_dict):
            potential_DMIs= []
            for dmi_id, dmi_matches in InterfaceHandling.protein_pairs_dict[protein_pair].dmi_matches_dict.items():
                for dmi_match in dmi_matches:
                    slim_protein= dmi_match.slim_protein
                    slim_match= f'{dmi_match.slim_match.start}-{dmi_match.slim_match.end}'
                    domain_protein= dmi_match.domain_protein
                    file.write('\t'.join((intx_ID, dmi_id, InterfaceHandling.SLiM_types_dict[dmi_id].name, InterfaceHandling.SLiM_types_dict[dmi_id].regex, dmi_match.slim_match.pattern, InterfaceHandling.SLiM_types_dict[dmi_id].probability, slim_protein, slim_match, domain_protein)))
                    for domain_match_list in dmi_match.domain_interface_match.domain_matches:
                        domain_id= domain_match_list[0].domain_id
                        start= [domain_match.start for domain_match in domain_match_list]
                        end= [domain_match.end for domain_match in domain_match_list]
                        domain_match= '|'.join([f"{i[0]}-{i[1]}" for i in zip(start, end)])
                        potential_DMIs.append(f'{slim_protein}:{dmi_id},{slim_match};{domain_protein}:{domain_id},{domain_match}')
                        print(f'{slim_protein}:{dmi_id},{slim_match};{domain_protein}:{domain_id},{domain_match}')
                        if domain_id in HMM_ref:
                            pfam_hmm= HMM_ref[domain_id]
                            if (domain_protein == protein_1) & (pfam_hmm in domains_1):
                                domains_1_overlap.add(f'{domain_id}:{pfam_hmm}')
                            elif (domain_protein == protein_2) & (pfam_hmm in domains_2):
                                domains_2_overlap.add(f'{domain_id}:{pfam_hmm}')
                        else:
                            if (domain_protein == protein_1) & (domain_id in domains_1):
                                domains_1_overlap.add(domain_id)
                            elif (domain_protein == protein_2) & (domain_id in domains_2):
                                domains_2_overlap.add(domain_id)
                        file.write('\t')
                        file.write('\t'.join((domain_id, domain_match)))
                    file.write('\n')
            df.loc[i, 'PotentialDMIMatch']= '||'.join([potential_DMI for potential_DMI in potential_DMIs])
            df.loc[i, 'domains_1_overlap']= ';'.join([domain for domain in domains_1_overlap])
            df.loc[i, 'domains_2_overlap']= ';'.join([domain for domain in domains_2_overlap])
            print('Potential DMI match annotated!')
    file.close()
    df.to_csv(interaction_file[:-4] + '_dmi_match_added.tsv', sep= '\t')
    print(f'New file saved as {interaction_file[:-4]}_dmi_match_added.tsv.')
    print(f'New file saved as {DMIMatch_outfile}.tsv.')

if __name__ == '__main__':

    prot_path= sys.argv[1]
    slim_type_file= sys.argv[2]
    dmi_type_file= sys.argv[3]
    smart_domain_types_file= sys.argv[4]
    pfam_domain_types_file= sys.argv[5]
    smart_domain_matches_json_file= sys.argv[6]
    pfam_domain_matches_json_file= sys.argv[7]
    features_path= sys.argv[8]
    interaction_file= sys.argv[9]
    DMIMatch_outfile= sys.argv[10]
    with open('/Users/chopyanlee/Coding/Python/DMI/domain_stuffs/smart_pfam_HMM_reference.txt', 'r') as f:
        HMM_ref= json.load(f)

    InterfaceHandling= DMIDB.InterfaceHandling(slim_type_file, dmi_type_file, smart_domain_types_file, pfam_domain_types_file, smart_domain_matches_json_file, pfam_domain_matches_json_file, features_path)
    InterfaceHandling.read_in_slim_types()
    InterfaceHandling.read_in_DMI_types()
    InterfaceHandling.read_in_domain_types()
    find_DMI_matches_in_interaction_file(prot_path, interaction_file, DMIMatch_outfile)

    # python3 find_DMI_match_in_interaction_file.py ~/Coding/Python/DMI/protein_sequences_and_features/human_protein_sequences ~/Coding/Python/DMI/elm_classes_20210222.tsv ~/Coding/Python/DMI/elm_interaction_domains_complete_20210222.tsv ~/Coding/Python/DMI/domain_stuffs/all_smart_domains_with_frequency.txt ~/Coding/Python/DMI/domain_stuffs/all_pfam_domains_with_frequency.txt ~/Coding/Python/DMI/domain_stuffs/interpro_9606_smart_matches_20210122.json ~/Coding/Python/DMI/domain_stuffs/interpro_9606_pfam_matches_20210122.json ~/Coding/Python/DMI/protein_sequences_and_features/human_protein_sequences_features interaction_data_1_20210309.tsv potential_DMI_in_CS_interaction_file_20210312
