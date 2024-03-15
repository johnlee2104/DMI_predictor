# This script uses the RRSFormation class from RRSFormation.py to form RRSv2. RRSv2 involves sampling a fixed number of DMI instances irrespective of DMI types from PPIs randomized using proteins from the PRS.
# Author: Chop Yan Lee

# import protein_interaction_interfaces
from DMIDB import *
import DMIDB
import RRSFormation
import sys, random, itertools

class RRSv2Formation(RRSFormation.RRSFormation):
    """
    Represents random reference set version 2

    Inherits from RRSFormation.RRSFormation
    """

    def __init__(self, RRS_version):
        """
        Instantiate RRS

        Args:
            RRS_version (str): a version name given to RRS, e.g. RRSv1_1_20210427 as RRS version 1, triplicate 1 and date of generating the RRS
        """
        super().__init__(RRS_version)

    def make_random_protein_pairs(self):
        """
        Generate random protein pairs using proteins from InterfaceHandling.proteins_dict and save the random pairs as ProteinPair instances. Additionally checks if the randomized protein pairs coincide with a known PPI, in which case the randomized protein pair is excluded.
        """
        print('Generating random pairs of proteins...')
        random_protein_pair= [tuple(sorted(pp)) for pp in itertools.combinations(list(InterfaceHandling.proteins_dict), 2)]
        print(f'{len(random_protein_pair)} combinations of protein pairs generated...')
        pairs_set= set(random_protein_pair)

        print('Removing known PPIs from randomly paired proteins...')
        filtered_pairs= pairs_set.difference(InterfaceHandling.known_PPIs)

        print('Creating ProteinPair instances...')
        for protein_pair in random.sample(list(filtered_pairs), 10000): # Sample only 10000 random protein pairs for DMI matching
            InterfaceHandling.protein_pairs_dict[protein_pair]= DMIDB.ProteinPair(protein_pair[0], protein_pair[1])
        print(f'{len(InterfaceHandling.protein_pairs_dict)} ProteinPair instances created.')

    def select_RRS_instances(self, number_instances):
        """
        Randomly sample a fixed number of instances irrespective of DMI types and append the sampled instances into self.RRS_instances

        Args:
            number_instances (int): Number of RRS instances to be sampled for each DMI type
        """
        sampled_dmi_matches= []

        print(f'Sampling DMI matches from {len(InterfaceHandling.protein_pairs_dict)} protein pairs...')
        while len(sampled_dmi_matches) < number_instances:
            for protpair in random.sample(list(InterfaceHandling.protein_pairs_dict), 1):
                protpair_inst= InterfaceHandling.protein_pairs_dict[protpair]
                if any(protpair_inst.dmi_matches_dict):
                    for slim_id in random.sample(list(protpair_inst.dmi_matches_dict), 1):
                        dmi_match_inst_list= protpair_inst.dmi_matches_dict[slim_id]
                        sampled_dmi_matches= sampled_dmi_matches + random.sample(dmi_match_inst_list, 1)
        print(f'Sampled {number_instances} DMI matches from {len(InterfaceHandling.protein_pairs_dict)} protein pairs.')

        self.RRS_instances= sampled_dmi_matches
        print(f'DMI match sampling completed. Number of DMI match sampled: {number_instances}.')

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
    number_instances= int(sys.argv[10])
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
        RRSv2= RRSv2Formation(RRS)
        InterfaceHandling.protein_pairs_dict= {}
        RRSv2.make_random_protein_pairs()
        InterfaceHandling.find_DMI_matches()
        RRSv2.select_RRS_instances(number_instances)
        RRSv2.write_out_RRS_instances(InterfaceHandling)

    # python3 RRSv2Formation.py ~/Coding/Python/DMI/protein_sequences_and_features/PRS_v3_RRS_v1_sequences PRS_hi_lit17_IntAct_known_PPIs_20210427.txt ../elm_classes_20210222.tsv ../elm_interaction_domains_complete_20210222.tsv ../domain_stuffs/all_smart_domains_with_frequency.txt ../domain_stuffs/all_pfam_domains_with_frequency.txt ../domain_stuffs/interpro_9606_smart_matches_20210122.json ../domain_stuffs/interpro_9606_pfam_matches_20210122.json ../protein_sequences_and_features/human_protein_sequences_features 1000 RRSv2_1_20210428 RRSv2_2_20210428 RRSv2_3_20210428
