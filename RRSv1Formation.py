import RRSFormation

class RRSv1Formation(RRSFormation.RRSFormation):
    def __init__(self, RRS_version):
        super().__init__(RRS_version)

    def make_random_protein_pairs(self): # Should I add an arg to specify how many RPP to make?
        pairs_set= set()

        for p in itertools.combinations(self.proteins_dict.keys(), 2):
            random_protein_pair= sorted(p)
            pairs_set.add(tuple(random_protein_pair))
        for protein_pair in pairs_set.difference(set(self.known_PPIs)):
            if protein_pair not in self.protein_pairs_dict:
                self.protein_pairs_dict[protein_pair] = ProteinPair(protein_pair[0], protein_pair[1])

    def select_RRS_instances(self, number_instances):
        for slim_id in self.DMI_types_dict.keys():
            dmi_matches= []
            for protpair, protpair_inst in self.protein_pairs_dict.items():
                if slim_id in protpair_inst.dmi_matches_dict.keys():
                    dmi_match_insts= protpair_inst.dmi_matches_dict[slim_id]
                    for random_inst in random.sample(dmi_match_insts, 1):
                        dmi_matches.append(random_inst)
            if len(dmi_matches) <= number_instances:
                self.RRS_instances = self.RRS_instances + dmi_matches
            else:
                selected_dmi_matches = random.sample(dmi_matches, number_instances)
                self.RRS_instances = self.RRS_instances + selected_dmi_matches

if __name__ == '__main__':

    domain_type_filename= sys.argv[1] # sys.argv[0] is the python script
    outputfile_name= sys.argv[2]

    DMI_DB= InterfaceHandling(files...)
    RRS_inst= RRSv1Formation('RRSv1')
    RRS_inst.DMIpredictor.read_in_domain_types(domain_type_filename)
    # Make sure to read in DMI types first, because we want to read in only domain types that are involved in DMI.
