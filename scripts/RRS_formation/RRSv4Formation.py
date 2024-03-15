# This script uses the RRSFormation class from RRSFormation.py to form RRSv4. RRSv4 involves sampling a fixed number of DMI instances irrespective of DMI types from PPIs randomized using all human proteins from SwissProt.
# Author: Chop Yan Lee

# For RRSv3 and 4, make sure to only sample from proteins with network information.
# import protein_interaction_interfaces
from DMIDB import *
import DMIDB
import RRSFormation
import sys, random, itertools

class RRSv4Formation(RRSFormation.RRSFormation):
    """
    Represents random reference set version 4

    Inherits from RRSFormation.RRSFormation
    """

    def __init__(self, RRS_version):
        """
        Instantiate RRS

        Args:
            RRS_version (str): a version name given to RRS, e.g. RRSv1_1_20210427 as RRS version 1, triplicate 1 and date of generating the RRS
        """
        super().__init__(RRS_version)

    def make_random_protein_pairs_with_network(self):
        """
        Generate random PPI using only proteins with network information
        """
        proteins_with_networks= set()

        file_names= [file_name for file_name in glob.glob(InterfaceHandling.network_path + '/*')]
        for file_name in file_names:
            prot_file= file_name.split('/')[-1]
            prot_id= prot_file.split('_')[0]
            proteins_with_networks.add(prot_id)

        print('Removing protein instances that do not have network information...')
        for prot in set(InterfaceHandling.proteins_dict).difference(proteins_with_networks): # remove prot instances without network info
            del InterfaceHandling.proteins_dict[prot]

        print('Generating protein combinations from proteins with network...')
        random_protein_pair= [tuple(sorted(pp)) for pp in itertools.combinations(list(InterfaceHandling.proteins_dict), 2)]
        print(f'{len(random_protein_pair)} combinations of protein pairs generated...')

        print('Creating ProteinPair instances...')
        for protein_pair in random.sample(random_protein_pair, 50000): # Sample only 50000 random protein pairs for DMI matching
            if protein_pair not in InterfaceHandling.known_PPIs:
                InterfaceHandling.protein_pairs_dict[protein_pair]= DMIDB.ProteinPair(protein_pair[0], protein_pair[1])
        print(f'{len(InterfaceHandling.protein_pairs_dict)} ProteinPair instances created.')

        sel_proteins= set()

        print('Removing protein instances that are not sampled as protein pair...')
        for protein_pair in InterfaceHandling.protein_pairs_dict:
            for protein in protein_pair:
                sel_proteins.add(protein)

        for prot in set(InterfaceHandling.proteins_dict).difference(sel_proteins):
            del InterfaceHandling.proteins_dict[prot]
        print(f'Total proteins sampled and paired as protein pair: {len(InterfaceHandling.proteins_dict)}.')
        
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
    network_path= sys.argv[10]
    number_instances= int(sys.argv[11])
    RRS_version= list(sys.argv[12:])

    InterfaceHandling= DMIDB.InterfaceHandling(prot_path, slim_type_file, dmi_type_file, smart_domain_types_file, pfam_domain_types_file, smart_domain_matches_json_file, pfam_domain_matches_json_file, features_path, network_path= network_path, PPI_file= PPI_file)
    InterfaceHandling.read_in_proteins()
    InterfaceHandling.read_in_known_PPIs()
    InterfaceHandling.read_in_slim_types()
    InterfaceHandling.read_in_DMI_types()
    InterfaceHandling.read_in_domain_types()
    InterfaceHandling.read_in_domain_matches()
    for RRS in RRS_version:
        RRSv4= RRSv4Formation(RRS)
        InterfaceHandling.protein_pairs_dict= {}
        RRSv4.make_random_protein_pairs_with_network()
        InterfaceHandling.create_slim_matches_all_proteins()
        InterfaceHandling.find_DMI_matches()
        RRSv4.select_RRS_instances(number_instances)
        RRSv4.write_out_RRS_instances(InterfaceHandling)

    # python3 RRSv4Formation.py ~/Coding/Python/DMI/protein_sequences_and_features/human_protein_sequences PRS_hi_lit17_IntAct_known_PPIs_20210427.txt ../elm_classes_20210222.tsv ../elm_interaction_domains_complete_20210222.tsv ../domain_stuffs/all_smart_domains_with_frequency.txt ../domain_stuffs/all_pfam_domains_with_frequency.txt ../domain_stuffs/interpro_9606_smart_matches_20210122.json ../domain_stuffs/interpro_9606_pfam_matches_20210122.json ../protein_sequences_and_features/human_protein_sequences_features ../protein_sequences_and_features/human_protein_sequences_features/Protein_networks_PRS_filtered 1000 RRSv4_1_20210428 RRSv4_2_20210428 RRSv4_3_20210428
