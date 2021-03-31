""" Takes a list of interaction pairs and find DMI matches between two interacting proteins. """

# This script generates two files, one is the original file annotated with potential DMI match, another is DMI match files, with similar format as PRS and RRS of DMI_DB. Note: start and end of SLiM and domain matches in annotated original file are annotated following positions of PDB structure, while those in DMI match files are annotated following canonical sequence from UniProt DB.
import sys, json, re
sys.path.append('/Users/chopyanlee/Coding/Python/DMI/protein_interaction_prediction/')
from DMIDB import *
import DMIDB
import pandas as pd

def find_DMI_matches_in_interaction_file(prot_path, interaction_file, DMIMatch_outfile):
    df= pd.read_csv(interaction_file, sep= '\t')
    file= open(DMIMatch_outfile + '.tsv', 'w')
    file.write('\t'.join(('intx_ID', 'Accession', 'Elm', 'Regex', 'Pattern', 'Probability', 'interactorElm', 'ElmMatch', 'interactorDomain', 'DomainID1', 'DomainMatch1', 'DomainID2', 'DomainMatch2')))
    file.write('\n')
    df['fraction_IUPlong_1']= float('NaN')
    df['fraction_IUPlong_2']= float('NaN')
    df['PotentialDMIMatch']= float('NaN')
    df['DMI_overlap_with_domains_1']= float('NaN')
    df['DMI_overlap_with_domains_2']= float('NaN')
    df['filter_by_interface']= float('NaN')
    df['domains_overlap_with_interface_1']= float('NaN')
    df['domains_overlap_with_interface_2']= float('NaN')
    for i, r in df.iterrows():
        protein_1= r['protein_1']
        protein_2= r['protein_2']
        intx_ID= r['intx_ID']
        seq_start_1= int(r['seq_start_1'])
        seq_end_1= int(r['seq_end_1'])
        seq_start_2= int(r['seq_start_2'])
        seq_end_2= int(r['seq_end_2'])
        IUPredLong_scores_1= []
        IUPredLong_scores_2= []
        try:
            with open(features_path + '/IUPred_long/' + protein_1 + '_iupredlong.txt', 'r') as f:
                lines= [line.strip() for line in f.readlines()]
            for line in lines[1:]:
                tabs= line.split('\t')
                IUPredLong_scores_1.append(float(tabs[2]))
            IUPredLong_scores_1= np.array(IUPredLong_scores_1)
            df.loc[i, 'fraction_IUPlong_1']= sum(IUPredLong_scores_1 >= 0.5)/len(IUPredLong_scores_1)
            print(f'fraction_IUPlong_1 of {protein_1} saved as {sum(IUPredLong_scores_1 >= 0.5)/len(IUPredLong_scores_1)}.')
        except:
            print(f'IUPredLong score of {protein_1} not found in {features_path}/IUPred_long/{protein_1}_iupredlong.txt.')
            print('Process terminated.')
            break
        try:
            with open(features_path + '/IUPred_long/' + protein_2 + '_iupredlong.txt', 'r') as f:
                lines= [line.strip() for line in f.readlines()]
            for line in lines[1:]:
                tabs= line.split('\t')
                IUPredLong_scores_2.append(float(tabs[2]))
            IUPredLong_scores_2= np.array(IUPredLong_scores_2)
            df.loc[i, 'fraction_IUPlong_2']= sum(IUPredLong_scores_2 >= 0.5)/len(IUPredLong_scores_2)
            print(f'fraction_IUPlong_2 of {protein_2} saved as {sum(IUPredLong_scores_2 >= 0.5)/len(IUPredLong_scores_2)}.')
        except:
            print(f'IUPredLong score of {protein_2} not found.')
            print('Process terminated.')
            break
        interface_residues_1= [int(re.sub('\D', '', x)) for x in r['res_iface_1'].split(';')]
        interface_residues_2= [int(re.sub('\D', '', x)) for x in r['res_iface_2'].split(';')]
        res_start_pdbID_1= 0
        res_end_pdbID_1= 0
        res_start_pdbID_2= 0
        res_end_pdbID_2= 0
        pos_diff_1= 0
        pos_diff_2= 0
        if (pd.notna(r['res_start_pdbID_1'])) & (pd.notna(r['res_end_pdbID_1'])):
            res_start_pdbID_1= int(r['res_start_pdbID_1'])
            res_end_pdbID_1= int(r['res_end_pdbID_1'])
        if (pd.notna(r['res_start_pdbID_2'])) & (pd.notna(r['res_end_pdbID_2'])):
            res_start_pdbID_2= int(r['res_start_pdbID_2'])
            res_end_pdbID_2= int(r['res_end_pdbID_2'])
        if (res_start_pdbID_1 - seq_start_1) != (res_end_pdbID_1 - seq_end_1):
            if (res_start_pdbID_1 != 0) & (res_end_pdbID_1 != 0):
                df.loc[i, 'filter_by_interface']= 'Start and end positions difference between pdbID and seq columns are not the same'
                continue
            else:
                if res_start_pdbID_1 != 0:
                    pos_diff_1= res_start_pdbID_1 - seq_start_1
                elif res_end_pdbID_1 != 0:
                    pos_diff_1= res_end_pdbID_1 - seq_end_1
        else:
            pos_diff_1= res_start_pdbID_1 - seq_start_1
        if (res_start_pdbID_2 - seq_start_2) != (res_end_pdbID_2 - seq_end_2):
            if (res_start_pdbID_2 != 0) & (res_end_pdbID_2 != 0):
                df.loc[i, 'filter_by_interface']= 'Start and end positions difference between pdbID and seq columns are not the same'
                continue
            else:
                if res_start_pdbID_2 != 0:
                    pos_diff_2= res_start_pdbID_2 - seq_start_2
                elif res_end_pdbID_2 != 0:
                    pos_diff_2= res_end_pdbID_2 - seq_end_2
        else:
            pos_diff_2= res_start_pdbID_2 - seq_start_2
        domains_1= set()
        domains_2= set()
        if pd.notna(r['domains_1']):
            if ';' in r['domains_1']:
                for domain in r['domains_1'].split(';'):
                    domains_1.add(domain)
            else:
                domains_1.add(r['domains_1'])
        if pd.notna(r['domains_2']):
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
            try:
                with open(prot_path + '/' + protein_id + '.txt', 'r') as f:
                    lines= [line.strip() for line in f.readlines()]
            except:
                print(f'{protein_id} sequence not found in {prot_path}')
                print('Process terminated.')
                break
            prot_inst= Protein(protein_id)
            InterfaceHandling.proteins_dict[protein_id]= prot_inst
            InterfaceHandling.proteins_dict[protein_id].sequence= lines[1]
        InterfaceHandling.create_slim_matches_all_proteins()
        for slim_id, slim_matches in InterfaceHandling.proteins_dict[protein_1].slim_matches_dict.items():
            filtered_slim_matches_1= []
            print(f'{protein_1} has {len(InterfaceHandling.proteins_dict[protein_1].slim_matches_dict[slim_id])} {slim_id} match before')
            for slim_match in slim_matches:
                if len(set(range(slim_match.start + pos_diff_1, slim_match.end + pos_diff_1 + 1)).intersection(set(interface_residues_1))) != 0:
                    filtered_slim_matches_1.append(slim_match)
            InterfaceHandling.proteins_dict[protein_1].slim_matches_dict[slim_id]= filtered_slim_matches_1
            print(f'{protein_1} has {len(InterfaceHandling.proteins_dict[protein_1].slim_matches_dict[slim_id])} {slim_id} match after')
        for slim_id, slim_matches in InterfaceHandling.proteins_dict[protein_2].slim_matches_dict.items():
            filtered_slim_matches_2= []
            print(f'{protein_2} has {len(InterfaceHandling.proteins_dict[protein_2].slim_matches_dict[slim_id])} {slim_id} match before')
            for slim_match in slim_matches:
                if len(set(range(slim_match.start + pos_diff_2, slim_match.end + pos_diff_2 + 1)).intersection(set(interface_residues_2))) != 0:
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
                if len(set(range(domain_match.start + pos_diff_1, domain_match.end + pos_diff_1)).intersection(set(interface_residues_1))) != 0:
                    filtered_domain_matches_1.append(domain_match)
            InterfaceHandling.proteins_dict[protein_1].domain_matches_dict[domain_id]= filtered_domain_matches_1
            print(f'{protein_1} has {len(InterfaceHandling.proteins_dict[protein_1].domain_matches_dict[domain_id])} {domain_id} match after')
        for domain_id, domain_matches in InterfaceHandling.proteins_dict[protein_2].domain_matches_dict.items():
            filtered_domain_matches_2= []
            print(f'{protein_2} has {len(InterfaceHandling.proteins_dict[protein_2].domain_matches_dict[domain_id])} {domain_id} match before')
            for domain_match in domain_matches:
                if len(set(range(domain_match.start + pos_diff_2, domain_match.end + pos_diff_2)).intersection(set(interface_residues_2))) != 0:
                    filtered_domain_matches_2.append(domain_match)
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
                    dmi_file_slim_match= f'{dmi_match.slim_match.start}-{dmi_match.slim_match.end}'
                    if slim_protein == protein_1:
                        slim_match= f'{dmi_match.slim_match.start + pos_diff_1}-{dmi_match.slim_match.end + pos_diff_1}'
                    elif slim_protein == protein_2:
                        slim_match= f'{dmi_match.slim_match.start + pos_diff_2}-{dmi_match.slim_match.end + pos_diff_2}'
                    domain_protein= dmi_match.domain_protein
                    file.write('\t'.join((intx_ID, dmi_id, InterfaceHandling.SLiM_types_dict[dmi_id].name, InterfaceHandling.SLiM_types_dict[dmi_id].regex, dmi_match.slim_match.pattern, InterfaceHandling.SLiM_types_dict[dmi_id].probability, slim_protein, dmi_file_slim_match, domain_protein)))
                    for domain_match_list in dmi_match.domain_interface_match.domain_matches:
                        domain_id= domain_match_list[0].domain_id
                        dmi_file_start= [domain_match.start for domain_match in domain_match_list]
                        dmi_file_end= [domain_match.end for domain_match in domain_match_list]
                        if domain_protein == protein_1:
                            start= [domain_match.start + pos_diff_1 for domain_match in domain_match_list]
                            end= [domain_match.end + pos_diff_1 for domain_match in domain_match_list]
                        elif domain_protein == protein_2:
                            start= [domain_match.start + pos_diff_2 for domain_match in domain_match_list]
                            end= [domain_match.end + pos_diff_2 for domain_match in domain_match_list]
                        domain_match= '|'.join([f'{i[0]}-{i[1]}' for i in zip(start, end)])
                        dmi_file_domain_match= '|'.join([f'{i[0]}-{i[1]}' for i in zip(dmi_file_start, dmi_file_end)])
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
                        file.write('\t'.join((domain_id, dmi_file_domain_match)))
                    file.write('\n')
            df.loc[i, 'PotentialDMIMatch']= '||'.join([potential_DMI for potential_DMI in potential_DMIs])
            df.loc[i, 'DMI_overlap_with_domains_1']= ';'.join([domain for domain in domains_1_overlap])
            df.loc[i, 'DMI_overlap_with_domains_2']= ';'.join([domain for domain in domains_2_overlap])
            print('||'.join([potential_DMI for potential_DMI in potential_DMIs]))
            if pd.isna(r['res_start_pdbID_1']) | pd.isna(r['res_end_pdbID_1']) | pd.isna(r['res_start_pdbID_2']) | pd.isna(r['res_end_pdbID_2']):
                df.loc[i, 'filter_by_interface']= 'DMI match found but sequence mapping might not be correct due to null value'
            else:
                df.loc[i, 'filter_by_interface']= 'Filtered by interface'
        domains_overlap_with_interface_1= []
        domains_overlap_with_interface_2= []
        for data in [smart_domain_matches, pfam_domain_matches]:
            for result in data['results']:
                if result['metadata']['accession']== protein_1:
                    for domain_match_id in result['entry_subset']:
                        for domain_match in domain_match_id['entry_protein_locations']:
                            for fragment in domain_match['fragments']:
                                if len(set(range(fragment['start'] + pos_diff_1, fragment['end'] + pos_diff_1)).intersection(set(interface_residues_1))) != 0:
                                    domains_overlap_with_interface_1.append(domain_match_id['accession'])
                elif result['metadata']['accession']== protein_2:
                    for domain_match_id in result['entry_subset']:
                        for domain_match in domain_match_id['entry_protein_locations']:
                            for fragment in domain_match['fragments']:
                                if len(set(range(fragment['start'] + pos_diff_2, fragment['end'] + pos_diff_2)).intersection(set(interface_residues_2))) != 0:
                                    domains_overlap_with_interface_2.append(domain_match_id['accession'])
        if any(domains_overlap_with_interface_1):
            df.loc[i, 'domains_overlap_with_interface_1']= ','.join([domain for domain in domains_overlap_with_interface_1])
            print('domains_overlap_with_interface_1 annotated.')
        if any(domains_overlap_with_interface_2):
            df.loc[i, 'domains_overlap_with_interface_2']= ','.join([domain for domain in domains_overlap_with_interface_2])
            print('domains_overlap_with_interface_2 annotated.')
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
    with open('/Users/chopyanlee/Coding/Python/DMI/domain_stuffs/interpro_9606_smart_matches_20210122.json', 'r') as f:
        smart_domain_matches= json.load(f)
    print('smart_json loaded')
    with open('/Users/chopyanlee/Coding/Python/DMI/domain_stuffs/interpro_9606_pfam_matches_20210122.json', 'r') as f:
        pfam_domain_matches= json.load(f)
    print('pfam_json loaded')

    InterfaceHandling= DMIDB.InterfaceHandling(slim_type_file, dmi_type_file, smart_domain_types_file, pfam_domain_types_file, smart_domain_matches_json_file, pfam_domain_matches_json_file, features_path)
    InterfaceHandling.read_in_slim_types()
    InterfaceHandling.read_in_DMI_types()
    InterfaceHandling.read_in_domain_types()
    find_DMI_matches_in_interaction_file(prot_path, interaction_file, DMIMatch_outfile)

    # python3 find_DMI_match_in_interaction_file.py ~/Coding/Python/DMI/protein_sequences_and_features/human_protein_sequences ~/Coding/Python/DMI/elm_classes_20210222.tsv ~/Coding/Python/DMI/elm_interaction_domains_complete_20210222.tsv ~/Coding/Python/DMI/domain_stuffs/all_smart_domains_with_frequency.txt ~/Coding/Python/DMI/domain_stuffs/all_pfam_domains_with_frequency.txt ~/Coding/Python/DMI/domain_stuffs/interpro_9606_smart_matches_20210122.json ~/Coding/Python/DMI/domain_stuffs/interpro_9606_pfam_matches_20210122.json ~/Coding/Python/DMI/protein_sequences_and_features/human_protein_sequences_features interaction_data_1_20210309.tsv potential_DMI_in_CS_interaction_file_20210331
