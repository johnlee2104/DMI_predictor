# import protein_interaction_interfaces
from DMIDB import *
import DMIDB
import RRSFormation
import sys, random, itertools

class RRSv1Formation(RRSFormation.RRSFormation):
    def __init__(self, RRS_version):
        super().__init__(RRS_version)

    def make_random_protein_pairs(self): # Should I add an arg to specify how many RPP to make?

        pairs_set= set()

        for p in itertools.combinations(InterfaceHandling.proteins_dict.keys(), 2):
            random_protein_pair= sorted(p)
            pairs_set.add(tuple(random_protein_pair))
        for protein_pair in pairs_set.difference(set(InterfaceHandling.known_PPIs)):
            if protein_pair not in InterfaceHandling.protein_pairs_dict:
                InterfaceHandling.protein_pairs_dict[protein_pair] = DMIDB.ProteinPair(protein_pair[0], protein_pair[1])

    def select_RRS_instances(self, number_instances):
        number_instances= int(number_instances)
        for slim_id in InterfaceHandling.DMI_types_dict.keys():
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

    # def write_out_RRS_instances(self): # take RRS_version as argument to construct the file name of the output file
    #     file= open(self.write_file, 'w')
    #     file.write('\t'.join((' ', 'Accession', 'Elm', 'Regex', 'Probability', 'interactorElm', 'ElmMatch', 'interactorDomain', 'DomainID1', 'DomainMatch1', 'DomainMatchEvalue1', 'DomainID2', 'DomainMatch2', 'DomainMatchEvalue2')))
    #     file.write('\n')
    #     for i , inst in enumerate(self.RRS_instances):
    #         file.write('\t'.join((str(i+1), inst.slim_match.slim_id, InterfaceHandling.SLiM_types_dict[inst.slim_match.slim_id].name, InterfaceHandling.SLiM_types_dict[inst.slim_match.slim_id].regex, InterfaceHandling.SLiM_types_dict[inst.slim_match.slim_id].probability, inst.slim_protein, '-'.join([str(inst.slim_match.start) , str(inst.slim_match.end)]), inst.domain_protein)))
    #         for domain_match_list in inst.domain_interface_match.domain_matches:
    #             domain_id= domain_match_list[0].domain_id
    #             start= [domain_match.start for domain_match in domain_match_list]
    #             end= [domain_match.end for domain_match in domain_match_list]
    #             evalue= [domain_match.evalue for domain_match in domain_match_list]
    #             match= '|'.join([f"{i[0]}-{i[1]}" for i in zip(start, end)])
    #             evalues= '|'.join([str(ev) for ev in evalue])
    #             file.write('\t')
    #             file.write('\t'.join((domain_id, match, evalues)))
    #         file.write('\n')
    #     file.close()
    #     print(f'{self.RRS_version} file saved as {self.write_file}.')

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
    RRS_version= sys.argv[10]
    number_instances= sys.argv[11]

    InterfaceHandling= DMIDB.InterfaceHandling(slim_type_file, dmi_type_file, smart_domain_types_file, pfam_domain_types_file, smart_domain_matches_json_file, pfam_domain_matches_json_file, features_path)
    # DMI_DB= DMIDB.InterfaceHandling(slim_type_file, dmi_type_file, smart_domain_types_file, pfam_domain_types_file, smart_domain_matches_json_file, pfam_domain_matches_json_file, features_path)
    RRS_v1= RRSv1Formation(RRS_version)
    InterfaceHandling.read_in_proteins(prot_path, canonical= True)
    InterfaceHandling.read_in_known_PPIs(PPI_file)
    RRS_v1.make_random_protein_pairs()
    InterfaceHandling.read_in_slim_types()
    InterfaceHandling.read_in_DMI_types()
    InterfaceHandling.read_in_domain_types()
    InterfaceHandling.read_in_domain_matches()
    InterfaceHandling.create_slim_matches_all_proteins()
    # InterfaceHandling.read_in_features_scores_all_proteins()
    # InterfaceHandling.calculate_average_features_scores_all_proteins()
    InterfaceHandling.find_DMI_matches()
    RRS_v1.select_RRS_instances(number_instances)
    RRS_v1.write_out_RRS_instances(InterfaceHandling)

    # python3 RRSv1Formation.py ~/Coding/Python/DMI/protein_sequences_and_features/PRS_v3_RRS_v1_sequences PRS_IntAct_union_known_PPIs.txt ../elm_classes_20210222.tsv ../elm_interaction_domains_complete_20210222.tsv ../domain_stuffs/all_smart_domains_with_frequency.txt ../domain_stuffs/all_pfam_domains_with_frequency.txt ../domain_stuffs/interpro_9606_smart_matches_20210122.json ../domain_stuffs/interpro_9606_pfam_matches_20210122.json ../protein_sequences_and_features/PRS_v3_RRS_v1_sequences_features RRSv1_1 5
