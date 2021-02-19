import sys
""" Comments... """
class Protein:
    def __init__(self, protein_id):
        self.protein_id= protein_id
        self.sequence= ''
        self.domain_matches_dict= {} # {'domain_id': [domain_match1, domain_match2, ...]}

class ProteinPair:
    def __init__(self, proteinA, proteinB):
        self.proteinA= proteinA
        self.proteinB= proteinB

class DomainType:
    def __init__(self, domain_id):
        self.domain_id= domain_id
        self.name= ''
        self.source= ''
        self.protein_frequency= None
        self.proteome_frequency= None

class DomainMatch:
    def __init__(self, domain_id, start, end):
        self.domain_id= domain_id
        self.start= start
        self.end= end
        self.evalue= None

class DomainInterface:
    def __init__(self):
        self.domain_dict= {}

class DomainInterfaceMatch:
    def __init__(self, domain_interface, domain_matches):
        self.domain_interface= domain_interface
        self.domain_matches= domain_matches

class InterfaceHandling:
    def __init__(self):
        self.proteins_dict= {}
        self.domain_types_dict= {}
        self.known_PPIs= []
        self.protein_pairs_dict= {}

    def read_in_proteins(self, in_file):
        with open(protein_id_seq_file,'r') as file:
            lines = [line.strip() for line in file.readlines()]
        for line in lines:
            if line[0] == '>':
                tab= line.split('|')
                protein_id= tab[1]
                prot_inst= Protein(protein_id)
                self.proteins_dict[protein_id]= prot_inst
            else:
                self.proteins_dict[protein_id].sequence += line

    def read_in_domain_types(self, in_file): # This one reads in all domain types, useful for DDI predictor
        with open(domain_types_file, 'r') as file:
            lines= [line.strip() for line in file.readlines()]
        for line in lines[2:]:
            tab= line.split('\t')
            name= tab[0]
            domain_id= tab[1]
            if len(tab) > 2:
                descr= tab[2]
            self.domain_types_dict[domain_id].name= name
            self.domain_types_dict[domain_id].descr= descr
            if domain_id[:2] == 'PF':
                self.domain_types_dict[domain_id].source= 'PFAM'
            elif domain_id[:2] == 'SM':
                self.domain_types_dict[domain_id].source= 'SMART'

    def read_in_known_PPIs(self, PPI_file):
        with open(PPI_file, 'r') as file:
            lines= [line.strip() for line in file.readlines()] # PRS saved as .tsv
        for line in lines:
            tab= line.split('\t')
            PPI_instance= sorted(list([tab[0], tab[1]])) # a sorted protein pair as list
            self.known_PPIs.append(tuple(PPI_instance)) # PPI pair saved as tuple

    def read_in_domain_matches(self): # This one reads in all domain matches, useful for DDI predictor
        with open(domain_matches_json_file) as f:
            data= json.load(f)
        for result in data['results']:
            protein_id= result['metadata']['accession']
            if protein_id in self.proteins_dict:
                self.proteins_dict[protein_id].name= result['metadata']['name']
                for domain_match_id in result['entry_subset']:
                    for domain_match in domain_match_id['entry_protein_locations']:
                        domain_id= domain_match['model']
                        score= domain_match['score']
                        for fragment in domain_match['fragments']:
                            start= fragment['start']
                            end= fragment['end']
                            domain_match_inst= DomainMatch(domain_id, start, end)
                            domain_match_inst.evalue= score
                            if domain_id not in self.proteins_dict[protein_id].domain_matches_dict.keys():
                                self.proteins_dict[protein_id].domain_matches_dict[domain_id]= []
                            self.proteins_dict[protein_id].domain_matches_dict[domain_id].append(domain_match_inst)