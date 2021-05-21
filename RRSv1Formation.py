# import protein_interaction_interfaces
from DMIDB import *
import DMIDB
import RRSFormation
import sys, random, itertools

class RRSv1Formation(RRSFormation.RRSFormation):
    def __init__(self, RRS_version):
        super().__init__(RRS_version)

    def make_random_protein_pairs(self):

        print('Generating random pairs of proteins...')
        random_protein_pair= [tuple(sorted(pp)) for pp in itertools.combinations(list(InterfaceHandling.proteins_dict), 2)]
        print(f'{len(random_protein_pair)} combinations of protein pairs generated...')
        pairs_set= set(random_protein_pair)

        print('Removing known PPIs from randomly paired proteins...')
        filtered_pairs= pairs_set.difference(InterfaceHandling.known_PPIs)

        print('Creating ProteinPair instances...')
        for protein_pair in filtered_pairs:
            InterfaceHandling.protein_pairs_dict[protein_pair]= DMIDB.ProteinPair(protein_pair[0], protein_pair[1])
        print(f'{len(InterfaceHandling.protein_pairs_dict)} ProteinPair instances created.')

    def select_RRS_instances(self, number_instances):
        number_instances= int(number_instances)
        for slim_id in InterfaceHandling.dmi_types_dict.keys():
            dmi_matches= []
            for protpair, protpair_inst in InterfaceHandling.protein_pairs_dict.items():
                if slim_id in protpair_inst.dmi_matches_dict.keys():
                    dmi_match_insts= protpair_inst.dmi_matches_dict[slim_id]
                    for random_inst in random.sample(dmi_match_insts, 1):
                        dmi_matches.append(random_inst)
            if len(dmi_matches) <= number_instances:
                self.RRS_instances = self.RRS_instances + dmi_matches
                print(f'{len(dmi_matches)} selected for {slim_id}.')
            else:
                selected_dmi_matches = random.sample(dmi_matches, number_instances)
                self.RRS_instances = self.RRS_instances + selected_dmi_matches
                print(f'{len(selected_dmi_matches)} selected for {slim_id}.')

if __name__ == '__main__':

    prot_path= sys.argv[1]
    PPI_file= sys.argv[2]
    slim_type_file= sys.argv[3]
    dmi_type_file= sys.argv[4]
    smart_domain_types_file= sys.argv[5]
    pfam_domain_types_file= sys.argv[6]
    smart_domain_matches_json_file= sys.argv[7]
    pfam_domain_matches_json_file= sys.argv[8]
    features_path= sys.argv[9]
    number_instances= sys.argv[10]
    RRS_version= list(sys.argv[11:])

    InterfaceHandling= DMIDB.InterfaceHandling(prot_path, slim_type_file, dmi_type_file, smart_domain_types_file, pfam_domain_types_file, smart_domain_matches_json_file, pfam_domain_matches_json_file, features_path, PPI_file= PPI_file)
    InterfaceHandling.read_in_proteins()
    InterfaceHandling.read_in_known_PPIs()
    InterfaceHandling.read_in_slim_types()
    InterfaceHandling.read_in_DMI_types()
    InterfaceHandling.read_in_domain_types()
    InterfaceHandling.read_in_domain_matches()
    InterfaceHandling.create_slim_matches_all_proteins()
    for RRS in RRS_version:
        RRS_v1= RRSv1Formation(RRS)
        RRS_v1.make_random_protein_pairs()
        InterfaceHandling.find_DMI_matches()
        RRS_v1.RRS_instances= []
        RRS_v1.select_RRS_instances(number_instances)
        RRS_v1.write_out_RRS_instances(InterfaceHandling)

    # python3 RRSv1Formation.py ~/Coding/Python/DMI/protein_sequences_and_features/PRS_v3_RRS_v1_sequences PRS_hi_lit17_IntAct_known_PPIs_20210427.txt ../elm_classes_20210222.tsv ../elm_interaction_domains_complete_20210222.tsv ../domain_stuffs/all_smart_domains_with_frequency.txt ../domain_stuffs/all_pfam_domains_with_frequency.txt ../domain_stuffs/interpro_9606_smart_matches_20210122.json ../domain_stuffs/interpro_9606_pfam_matches_20210122.json ../protein_sequences_and_features/human_protein_sequences_features 5 RRSv1_1_20210427 RRSv1_2_20210427 RRSv1_3_20210427
