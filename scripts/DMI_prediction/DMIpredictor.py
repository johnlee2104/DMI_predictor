# This script executes DMI prediction using classes and functions from DMIDB.py.
# Author: Chop Yan Lee

# Workflow:
## - read in DMI file and SLiM types file (ELM class with regex)
## - user input either a pair of interacting proteins or a list of interacting protein pair in .txt file. (optparse)
## - set up a protein object for each of the input protein and a ProteinPair object for each of the input proteinpair
## - read in the sequences and sequence features like IUPred, etc and domain matches for all the protein objects
## - for each ProteinPair object, if there is a SLiM-binding domain, run regex search of the SLiM on the other interacting protein. Do DMI matching.
## - calculate the SLiM feature for DMI matches that are detected by specifying the domain type option, missing features as np.nan so that they can be imputed
## - impute missing features if needed, then do prediction on the feature. The predicted probability should be saved in DMIMatch.score
## - output result in .tsv format similar to PRS and RRS. Add intx_id (protA_protB, sorted), which protein in the real network has the dommain match, degree of slim protein in the real network, comments (flag cases with domain count > 1, and those with two domains from two proteins), imputed feature
## - output should also be saved into MySQL

from DMIDB import *
import DMIDB
from optparse import OptionParser
import joblib, sys

parser= OptionParser()
parser.add_option('-p', '--pair', dest= 'proteinpair', help= 'Protein pair mode')
parser.add_option('-l', '--list', dest= 'proteinpairlist', help= 'Protein pair list mode')

(options, args)= parser.parse_args()

proteinpair= options.proteinpair
proteinpairlist= options.proteinpairlist

features= ['Probability', 'IUPredShort', 'Anchor', 'DomainOverlap', 'qfo_RLC', 'qfo_RLCvar', 'vertebrates_RLC', 'vertebrates_RLCvar', 'mammalia_RLC', 'mammalia_RLCvar', 'metazoa_RLC', 'metazoa_RLCvar', 'DomainEnrichment_pvalue', 'DomainEnrichment_zscore', 'DomainFreqbyProtein', 'DomainFreqinProteome']

def get_protein_pair():
    """
    Depending on the option, read the protein pair(s) into InterfaceHandling.protein_pairs_dict as DMIDB.ProteinPair instances

    Returns:
        input_proteins (set): UniProt IDs of individual proteins stored in a set
    """
    input_proteins= set()
    if proteinpair != None:
        protein_pair= tuple(sorted(proteinpair.split(',')))
        InterfaceHandling.protein_pairs_dict[protein_pair]= DMIDB.ProteinPair(protein_pair[0], protein_pair[1])
        for protein in protein_pair:
            input_proteins.add(protein)
    else:
        with open(proteinpairlist, 'r') as f:
            lines= [line.strip() for line in f.readlines()]
        for line in lines:
            protein_pair= tuple(sorted(line.split('\t')))
            InterfaceHandling.protein_pairs_dict[protein_pair]= DMIDB.ProteinPair(protein_pair[0], protein_pair[1])
            for protein in protein_pair:
                input_proteins.add(protein)
    return input_proteins

def create_slim_match_in_protein_pair(): # use the domain match in the proteins of a protein pair to restrict create_slim_matches to only relevant dmi types
    """
    Perform DMI matching by first checking the domain matches in proteinA, following by calling the method Protein.create_slim_matches in proteinB to find matches of SLiM types that can bind to the domain matches in proteinA. The same process is repeated by switching proteinA and B to find all DMI matches. As opposed to InterfaceHandling.create_slim_matches_all_proteins, this function restricts DMI search to only the SLiM types that can bind to any domains of the interacting protein.
    """
    for protpair in InterfaceHandling.protein_pairs_dict:
        protein_pairs= [(protpair[0], protpair[1]), (protpair[1], protpair[0])]
        for protein_pair in protein_pairs:
            prot_inst= InterfaceHandling.proteins_dict[protein_pair[0]]
            domain_matches_dict= prot_inst.domain_matches_dict
            sel_dmi_type_inst= []
            for dmi_type, dmi_type_inst in InterfaceHandling.dmi_types_dict.items():
                domain_interface= dmi_type_inst.domain_interfaces
                if len(domain_interface) == 1: # Cases of 1 domain or 2 domains in one protein
                    domains= set(dmi_type_inst.domain_interfaces[0].domain_dict.keys())
                    if domains.intersection(set(domain_matches_dict.keys())) == domains:
                        sel_dmi_type_inst.append(tuple([dmi_type_inst, domains]))
                else: # cases of len(cognate_domains) == 2 -> 2 domains in 2 proteins
                    domains= set()
                    for domain_int in domain_interface:
                        for domain_id in domain_int.domain_dict.keys():
                            domains.add(domain_id)
                        sel_domain= domains.intersection(set(domain_matches_dict.keys()))
                        if any(sel_domain):
                            sel_dmi_type_inst.append(tuple([dmi_type_inst, sel_domain]))
            for dmi_type_inst, domain_type_list in sel_dmi_type_inst:
                slim_type_inst= InterfaceHandling.slim_types_dict[dmi_type_inst.slim_id]
                InterfaceHandling.proteins_dict[protein_pair[1]].create_slim_matches(dmi_type_inst, slim_type_inst)
                if dmi_type_inst.slim_id in InterfaceHandling.proteins_dict[protein_pair[1]].slim_matches_dict:
                    for slim_match_inst in InterfaceHandling.proteins_dict[protein_pair[1]].slim_matches_dict[dmi_type_inst.slim_id]:
                        slim_match_inst.get_slim_match_features(domain_type_list)

def predict_DMI_match(model,imputer):
    """
    Iterate over InterfaceHandling.protein_pairs_dict and scores the DMI matches saved in the ProteinPair instances using a trained random forest model. The scores are saved in DMIMatch.score.

    Args:
        model (sklearn.ensemble.RandomForestClassifier): A trained random forest model
        imputer (sklearn.impute.SimpleImputer): A fitted imputer
    """
    for prot_pair, prot_pair_inst in InterfaceHandling.protein_pairs_dict.items():
        for slim_id, dmi_match_inst_list in prot_pair_inst.dmi_matches_dict.items():
            for dmi_match_inst in dmi_match_inst_list:
                feature_array= []
                feature_array.append(InterfaceHandling.slim_types_dict[slim_id].probability)
                slim_match_inst= dmi_match_inst.slim_match
                features_dict= slim_match_inst.__dict__
                for feature in list(features_dict)[6:19]:
                    if features_dict[feature] == None:
                        features_dict[feature]= ''
                        feature_array.append(np.nan)
                    else:
                        feature_array.append(features_dict[feature])
                domainfreqsbyprotein= []
                domainfreqsinproteome= []
                for domain_match_list in dmi_match_inst.domain_interface_match.domain_matches:
                    domain_id= domain_match_list[0].domain_id
                    domainfreqsbyprotein.append(float(InterfaceHandling.domain_types_dict[domain_id].DomainFreqbyProtein))
                    domainfreqsinproteome.append(float(InterfaceHandling.domain_types_dict[domain_id].DomainFreqinProteome))
                if len(domainfreqsbyprotein) > 1:
                    feature_array.append(np.mean(domainfreqsbyprotein))
                    feature_array.append(np.mean(domainfreqsinproteome))
                else:
                    feature_array= feature_array + domainfreqsbyprotein
                    feature_array= feature_array + domainfreqsinproteome
                feature_array= np.array(feature_array).reshape(1,-1)
                feature_array= feature_array.astype(float)
                transformed_feature_array= imputer.transform(feature_array)
                missing_indicator= np.isnan(feature_array)
                dmi_match_inst.score= model.predict_proba(transformed_feature_array)[:,1][0]
                dmi_match_inst.missing_feature= ','.join((np.array(features)[missing_indicator.ravel()]))

def write_out_DMI_match():
    """
    Write out DMI matches, as well as their features and scores, in a .tsv file. Provide additional information, such as additional binding domain from another partner is required for DMI to form, is also saved in the Notes column. 
    """
    if proteinpair != None:
        output_name= '_'.join([i for i in sorted(proteinpair.split(','))])
    else:
        output_name, _= proteinpairlist.split('.')
    with open(f'{output_name}_DMI_prediction.tsv', 'w') as f:
        f.write('\t'.join(('intx_ID', 'Accession', 'Elm', 'Regex', 'Pattern', 'Probability', 'SLiMProtein', 'SLiMMatch', 'IUPredShort', 'Anchor', 'DomainOverlap', 'qfo_RLC', 'qfo_RLCvar', 'vertebrates_RLC', 'vertebrates_RLCvar', 'mammalia_RLC', 'mammalia_RLCvar', 'metazoa_RLC', 'metazoa_RLCvar', 'DomainEnrichment_pvalue', 'DomainEnrichment_zscore', 'Partner_with_domain', 'DomainProtein', 'DomainID1', 'DomainName1','DomainMatch1', 'DomainMatchFound1', 'DomainMatchRequired1', 'DomainMatchEvalue1', 'DomainFreqbyProtein1', 'DomainFreqinProteome1', 'DomainID2', 'DomainName2', 'DomainMatch2', 'DomainMatchFound2', 'DomainMatchRequired2', 'DomainMatchEvalue2', 'DomainFreqbyProtein2', 'DomainFreqinProteome2', 'DMIMatchScore', 'Notes')))
        f.write('\n')
        for prot_pair, prot_pair_inst in InterfaceHandling.protein_pairs_dict.items():
            intx_ID= '_'.join([id for id in prot_pair])
            for slim_id, dmi_match_inst_list in prot_pair_inst.dmi_matches_dict.items():
                for dmi_match_inst in dmi_match_inst_list:
                    f.write(f'{intx_ID}\t')
                    f.write('\t'.join((dmi_match_inst.slim_match.slim_type_inst.slim_id, dmi_match_inst.slim_match.slim_type_inst.name, dmi_match_inst.slim_match.slim_type_inst.regex, dmi_match_inst.slim_match.pattern, dmi_match_inst.slim_match.slim_type_inst.probability, dmi_match_inst.slim_protein, '-'.join([str(dmi_match_inst.slim_match.start), str(dmi_match_inst.slim_match.end)]), str(dmi_match_inst.slim_match.IUPredShort), str(dmi_match_inst.slim_match.Anchor), str(dmi_match_inst.slim_match.DomainOverlap), str(dmi_match_inst.slim_match.qfo_RLC), str(dmi_match_inst.slim_match.qfo_RLCvar), str(dmi_match_inst.slim_match.vertebrates_RLC), str(dmi_match_inst.slim_match.vertebrates_RLCvar), str(dmi_match_inst.slim_match.mammalia_RLC), str(dmi_match_inst.slim_match.mammalia_RLCvar), str(dmi_match_inst.slim_match.metazoa_RLC), str(dmi_match_inst.slim_match.metazoa_RLCvar), str(dmi_match_inst.slim_match.DomainEnrichment_pvalue), str(dmi_match_inst.slim_match.DomainEnrichment_zscore), ','.join([p.protein_id for p in dmi_match_inst.slim_match.partners_with_domain_in_real_network]), dmi_match_inst.domain_protein)))
                    for domain_match_list in dmi_match_inst.domain_interface_match.domain_matches:
                        domain_id= domain_match_list[0].domain_id
                        domain_name= InterfaceHandling.domain_types_dict[domain_id].name
                        start= [domain_match.start for domain_match in domain_match_list]
                        end= [domain_match.end for domain_match in domain_match_list]
                        domain_count_found= len(domain_match_list)
                        for intf in InterfaceHandling.dmi_types_dict[slim_id].domain_interfaces:
                            if domain_id in intf.domain_dict.keys():
                                domain_count_required= intf.domain_dict[domain_id]
                        evalue= [domain_match.evalue for domain_match in domain_match_list]
                        match= '|'.join([f"{i[0]}-{i[1]}" for i in zip(start, end)])
                        evalues= '|'.join([str(ev) for ev in evalue])
                        domainfreqbyprotein= InterfaceHandling.domain_types_dict[domain_id].DomainFreqbyProtein
                        domainfreqinproteome= InterfaceHandling.domain_types_dict[domain_id].DomainFreqinProteome
                        f.write('\t')
                        f.write('\t'.join((domain_id, domain_name, match, str(domain_count_found), str(domain_count_required), evalues, domainfreqbyprotein, domainfreqinproteome)))
                    if len(dmi_match_inst.domain_interface_match.domain_matches) == 1:
                        f.write(f'\t\t\t\t\t\t\t\t')
                    f.write(f'\t{dmi_match_inst.score}')
                    if len(InterfaceHandling.dmi_types_dict[slim_id].domain_interfaces) > 1: # tripartite complexes
                        domain_required= []
                        for intf in InterfaceHandling.dmi_types_dict[slim_id].domain_interfaces:
                            domain_required= [(d, count) for d, count in intf.domain_dict.items() if d != domain_id]
                        comment= f'Binding of motif requires an additional partner with {domain_required[0][1]} match(es) of {domain_required[0][0]} ({InterfaceHandling.domain_types_dict[domain_required[0][0]].name}).'
                        f.write(f'\t{comment}')
                    # f.write(f'\t{dmi_match_inst.domain_interface_match.domain_interface.domain_dict}, {[intf.domain_dict for intf in InterfaceHandling.dmi_types_dict[slim_id].domain_interfaces]}')
                    f.write('\n')

if __name__ == '__main__':
    prot_path= sys.argv[1]
    slim_type_file= sys.argv[2]
    dmi_type_file= sys.argv[3]
    smart_domain_types_file= sys.argv[4]
    pfam_domain_types_file= sys.argv[5]
    smart_domain_matches_json_file= sys.argv[6]
    pfam_domain_matches_json_file= sys.argv[7]
    features_path= sys.argv[8]
    RF_path= sys.argv[9]
    RF= joblib.load(RF_path)
    imputer_path= sys.argv[10]
    imputer= joblib.load(imputer_path)

    InterfaceHandling= DMIDB.InterfaceHandling(prot_path, slim_type_file, dmi_type_file, smart_domain_types_file, pfam_domain_types_file, smart_domain_matches_json_file, pfam_domain_matches_json_file, features_path)
    InterfaceHandling.read_in_proteins()
    InterfaceHandling.read_in_slim_types()
    InterfaceHandling.read_in_DMI_types()
    InterfaceHandling.read_in_domain_types()
    InterfaceHandling.read_in_domain_matches()
    input_proteins= get_protein_pair()
    InterfaceHandling.read_in_networks(prot_set= input_proteins)
    for protein in input_proteins:
        InterfaceHandling.proteins_dict[protein].read_in_features(features_path)
    create_slim_match_in_protein_pair()
    InterfaceHandling.find_DMI_matches()
    predict_DMI_match(model=RF,imputer=imputer)
    write_out_DMI_match()

# obsolete command
# python3 DMIpredictor.py  ~/Coding/Python/DMI/protein_sequences_and_features/human_protein_sequences ../elm_classes_20210222.tsv ../elm_interaction_domains_complete_20210222.tsv ../domain_stuffs/all_smart_domains_with_frequency.txt ../domain_stuffs/all_pfam_domains_with_frequency.txt ../domain_stuffs/interpro_9606_smart_matches_20210122.json ../domain_stuffs/interpro_9606_pfam_matches_20210122.json ../protein_sequences_and_features/human_protein_sequences_features ~/Coding/Python/DMI/final_RF_model_with_RRSv4_3.joblib ~/Coding/Python/DMI/final_median_imputer_with_RRSv4_3.joblib -l test_protein_pair.txt

# use this as an example!
# python3 ~/Coding/Python/DMI/DMI_predictor/DMIpredictor.py  ~/Coding/Python/DMI/protein_sequences_and_features/human_protein_sequences ~/Coding/Python/DMI/DMI_types_annotation/20220311_elm_classes.tsv ~/Coding/Python/DMI/DMI_types_annotation/20220311_elm_interaction_domains_complete.tsv ~/Coding/Python/DMI/domain_stuffs/all_smart_domains_with_frequency.txt ~/Coding/Python/DMI/domain_stuffs/all_pfam_domains_with_frequency.txt ~/Coding/Python/DMI/domain_stuffs/interpro_9606_smart_matches_20210122.json ~/Coding/Python/DMI/domain_stuffs/interpro_9606_pfam_matches_20210122.json ~/Coding/Python/DMI/protein_sequences_and_features/human_protein_sequences_features ~/Coding/Python/DMI/ML_stuffs/final_RF_model_with_RRSv4_3.joblib ~/Coding/Python/DMI/ML_stuffs/final_median_imputer_with_RRSv4_3.joblib -l test_protein_pair.tx