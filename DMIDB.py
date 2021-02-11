import protein_interaction_interfaces

class Protein(protein_interaction_interfaces.Protein):
    def __init__(self, protein_id):
        super().__init__(protein_id)
        self.slim_matches_dict= {} # {'slim_id': [slim_match1, slim_match2, ...]}
        self.IUPLong_scores= []
        self.IUPShort_scores= []
        self.Anchor_scores= []

    def create_slim_matches(self):
        for slim_id, slim_type_inst in InterfaceHandling.SLiM_types_dict.items():
            slim_start= [match.start() + 1 for match in re.finditer('(?=(' + slim_type_inst.regex + '))', self.sequence)]
            match_pattern= [match.group(1) for match in re.finditer('(?=(' + slim_type_inst.regex + '))', self.sequence)]
            match_results= list(zip(slim_start, match_pattern))
            if len(match_results) >0:
                self.slim_matches_dict[slim_id]= []
                for match in match_results:
                    slim_match_inst= SLiMMatch(slim_id, match[0], match[0] + len(match[1]) - 1)
                    self.slim_matches_dict[slim_id].append(slim_match_inst)

    def read_in_iupred_score(self):
        pass

    def calculate_average_iupred_anchor(self): # Dont pre-calculate this for now, only run this function when we want to make a prediction
        for slim_id, slim_match in self.slim_matches_dict.items():
            for slim_match_inst in slim_match:
                start= int(slim_match_inst.start)
                end= int(slim_match_inst.end)
                slim_match_inst.IUPLong= float(sum(self.IUPLong_scores[start-1:end])/(end - start + 1))
                slim_match_inst.IUPShort= float(sum(self.IUPShort_scores[start-1:end])/(end - start + 1))
                slim_match_inst.Anchor= float(sum(self.Anchor_scores[start-1:end])/(end - start + 1))

class ProteinPair(protein_interaction_interfaces.ProteinPair):
    def __init__(self, proteinA, proteinB):
        super().__init__(proteinA, proteinB)
        self.dmi_matches_dict= {} # {'slim_id':[dmi_match1, dmi_match2, ...]}

class SLiMType:
    def __init__(self, slim_id):
        self.slim_id= slim_id
        self.name= ''
        self.regex= ''
        self.probability= ''

class SLiMMatch:
    def __init__(self, slim_id, start, end):
        self.slim_id= slim_id
        self.start= start
        self.end= end
        self.IUPLong= None
        self.IUPredShort= None
        self.Anchor= None

class DMIType:
    def __init__(self, slim_id):
        self.slim_id= slim_id
        self.domain_interfaces= [] # [{'domain_id1': 1, 'domain_id2': 10}] / [{'domain_id3':1}, {'domain_id4': 4}]

class DMIMatch:
    def __init__(self,slim_protein,domain_protein,slim_match,domain_interface_match):
        self.slim_protein = slim_protein
        self.domain_protein = domain_protein
        self.slim_match = slim_match
        self.domain_interface_match = domain_interface_match
        self.score = None

### Rework this part
class InterfaceHandling(protein_interaction_interfaces.InterfaceHandling):
    def __init__(self, slim_type_file, dmi_type_file, domain_type_file, domain_matches_json_file):
        self.SLiM_types_dict= {}
        self.DMI_types_dict= {}
        self.slim_type_file= slim_type_file
        self.dmi_type_file= dmi_type_file
        self.domain_type_file= domain_type_file
        self.domain_matches_json_file= domain_matches_json_file

    def read_in_slim_types(self):

        """ For some reason, the regex that starts with [] python automatically adds a "" around the string, so one needs to modify these manually """

        with open(self.slim_type_file, 'r') as file:
            lines = [line.strip() for line in file.readlines()]
        for line in lines[1:]:
            tab= line.split('\t')
            slim_type_inst= SLiMType(tab[0])
            slim_type_inst.name= tab[1]
            slim_type_inst.regex= tab[4].replace('"', '')
            slim_type_inst.probability= float(tab[5])
            self.SLiM_types_dict[tab[0]] = slim_type_inst

    def read_in_DMI_types(self):
        with open(self.dmi_type_file, 'r') as file:
            lines= [line.strip() for line in file.readlines()]
        for line in lines[1:]:
            tab= line.split('\t')
            if tab[6] == '1': # Default use 1 == 1
                slim_id= tab[0]
                domain_id= tab[2]
                domain_id2= ''
                domain_count= int(tab[5])
                DMI_type_inst= DMIType(slim_id)
                if slim_id not in self.DMI_types_dict.keys():
                    self.DMI_types_dict[slim_id] = DMI_type_inst
                if len(tab) < 11 or (len(tab) > 10 and tab[11] == '') :
                    domain_interface_inst= DomainInterface()
                    domain_interface_inst.domain_dict[domain_id]= domain_count
                elif (len(tab) > 7) & (tab[11] == '1'):
                    domain_id2= tab[7]
                    domain_count2= tab[10]
                    domain_interface_inst= DomainInterface()
                    domain_interface_inst.domain_dict[domain_id] = domain_count
                    domain_interface_inst.domain_dict[domain_id2] = domain_count2
                    if domain_id2 not in self.domain_types_dict:
                        self.domain_types_dict[domain_id2] = DomainType(domain_id2)
                    self.domain_types_dict[domain_id2].dmi_types.append(self.DMI_types_dict[slim_id])
                elif (len(tab) > 7) & (tab[11] == '0'):
                    continue
                self.DMI_types_dict[slim_id].domain_interfaces.append(domain_interface_inst)
                if domain_id not in self.domain_types_dict:
                    self.domain_types_dict[domain_id] = DomainType(domain_id)
                self.domain_types_dict[domain_id].dmi_types.append(self.DMI_types_dict[slim_id])

    def read_in_domain_type(self): # this one reads in only domains involved in DMI
        with open(self.domain_types_file, 'r') as file:
            lines= [line.strip() for line in file.readlines()]
        for line in lines[2:]:
            tab= line.split('\t')
            name= tab[0]
            domain_id= tab[1]
            if len(tab) > 2:
                descr= tab[2]
            if domain_id not in self.domain_types_dict:
                continue
            self.domain_types_dict[domain_id].name= name
            self.domain_types_dict[domain_id].descr= descr
            if domain_id[:2] == 'PF':
                self.domain_types_dict[domain_id].source= 'PFAM'
            elif domain_id[:2] == 'SM':
                self.domain_types_dict[domain_id].source= 'SMART'

    def read_in_domain_matches(self): # this one reads in only domain matches involved in DMI
        if len(self.domain_types_dict) == 0:
            self.read_in_domain_types()
            print('Warning: No DMI types have been read in before')
            print('Reading in domain types now...')
        with open(self.domain_matches_json_file) as f:
            data= json.load(f)
        for result in data['results']:
            protein_id= result['metadata']['accession']
            if protein_id in self.proteins_dict:
                self.proteins_dict[protein_id].name= result['metadata']['name']
                for domain_match_id in result['entry_subset']:
                    for domain_match in domain_match_id['entry_protein_locations']:
                        domain_id= domain_match['model']
                        score= domain_match['score']
                        if domain_id in self.domain_types_dict:
                            for fragment in domain_match['fragments']:
                                start= fragment['start']
                                end= fragment['end']
                                domain_match_inst= DomainMatch(domain_id, start, end)
                                domain_match_inst.evalue= score
                                if domain_id not in self.proteins_dict[protein_id].domain_matches_dict.keys():
                                    self.proteins_dict[protein_id].domain_matches_dict[domain_id]= []
                                self.proteins_dict[protein_id].domain_matches_dict[domain_id].append(domain_match_inst)

    def create_slim_matches_all_proteins(self):
        for prot_id, prot_inst in self.proteins_dict.items():
            prot_inst.create_slim_matches()

# read_in_iupred_score and calculate_average_iupred can go into Protein class too that takes one protein and run the function using the protein
# and in InterfaceHandling, make another function that runs the Protein.functions globally.
    def read_in_iupred_score_all_proteins(self):
        pass

    def calculate_average_iupred_anchor_all_proteins(self): # Dont pre-calculate this for now, only run this function when we want to make a prediction
        for prot_id, prot_inst in self.proteins_dict.items():
            prot_inst.calculate_average_iupred_anchor()

    def find_DMI_matches(self):
        for protpair, protpair_inst in self.protein_pairs_dict.items():
            protein_pairs= [(protpair[0], protpair[1]), (protpair[1], protpair[0])]
            for protein_pair in protein_pairs:
                prot1_domain_matches_dict= self.proteins_dict[protein_pair[0]].domain_matches_dict
                prot2_slim_matches_dict= self.proteins_dict[protein_pair[1]].slim_matches_dict
                unique_slim_ids = set()
                for domain_id in prot1_domain_matches_dict.keys():
                    slim_id_list= [dmi_type.slim_id for dmi_type in self.domain_types_dict[domain_id].dmi_types]
                    slim_id_match_list= set(list(filter(lambda slim_id: slim_id in prot2_slim_matches_dict.keys(), slim_id_list))) # Check if protB has the slim_id of a domain in protA
                    unique_slim_ids = unique_slim_ids.union(slim_id_match_list)
                for slim_id_match in unique_slim_ids:
                    domain_interface_list= self.DMI_types_dict[slim_id_match].domain_interfaces
                    for domain_interface in domain_interface_list:
                        cognate_domains= list(domain_interface.domain_dict.keys())
                        if set(cognate_domains).intersection(set(prot1_domain_matches_dict.keys())) == set(cognate_domains):
                            if len(cognate_domains) == 1:
                                domain_interface_match_inst= DomainInterfaceMatch(slim_id_match, domain_interface, [prot1_domain_matches_dict[cognate_domains[0]]])
                                domain_interface_match_inst.domain_matches_count.append(len(prot1_domain_matches_dict[cognate_domains[0]]))
                                print(f'{protein_pair[0]} has {cognate_domains[0]} and {protein_pair[1]} has {slim_id_match}')
                                for slim_match in prot2_slim_matches_dict[slim_id_match]:
                                    DMIMatch_inst= DMIMatch(protein_pair[1], protein_pair[0], slim_match, domain_interface_match_inst)
                                    if slim_id_match not in protpair_inst.dmi_matches_dict.keys():
                                        protpair_inst.dmi_matches_dict[slim_id_match] = []
                                    protpair_inst.dmi_matches_dict[slim_id_match].append(DMIMatch_inst)
                            elif len(cognate_domains) == 2:
                                domain_match1= prot1_domain_matches_dict[cognate_domains[0]] # a list
                                domain_match2= prot1_domain_matches_dict[cognate_domains[1]] # a list
                                domain_interface_match_inst= DomainInterfaceMatch(slim_id_match, domain_interface, [domain_match1, domain_match2])
                                domain_interface_match_inst.domain_matches_count = [len(domain_match1), len(domain_match2)]
                                print(f'{protein_pair[0]} has {cognate_domains[0]}, {cognate_domains[1]} and {protein_pair[1]} has {slim_id_match}')
                                for slim_match in prot2_slim_matches_dict[slim_id_match]:
                                    DMIMatch_inst= DMIMatch(protein_pair[1], protein_pair[0], slim_match, domain_interface_match_inst)
                                    if slim_id_match not in protpair_inst.dmi_matches_dict.keys():
                                        protpair_inst.dmi_matches_dict[slim_id_match] = []
                                    protpair_inst.dmi_matches_dict[slim_id_match].append(DMIMatch_inst)
