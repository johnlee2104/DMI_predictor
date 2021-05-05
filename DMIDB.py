import protein_interaction_interfaces
from protein_interaction_interfaces import *
import re, json, requests, sys, glob
import scipy.special as sc
import numpy as np

dummy_value= 88888

class Protein(protein_interaction_interfaces.Protein):
    def __init__(self, protein_id):
        super().__init__(protein_id)
        self.slim_matches_dict= {} # {'slim_id': [slim_match1, slim_match2, ...]}
        self.IUPredLong_scores= []
        self.IUPredShort_scores= []
        self.Anchor_scores= []
        self.DomainOverlap_scores= []
        self.qfo_RLC_scores= {} # follows the cons score format from ND, {'pos1':'score1', 'pos2':'score2'}
        self.vertebrates_RLC_scores= {}
        self.mammalia_RLC_scores= {}
        self.metazoa_RLC_scores= {}
        self.networks= {} # saved as {0:[prot_inst1,prot_inst2],1:[prot_inst3,prot_inst4]} with 0 being real network and 1-1000 being random network
        self.network_degree= None

    def create_slim_matches(self, dmi_type_inst, slim_type_inst):
        slim_start= [match.start() for match in re.finditer('(?=(' + slim_type_inst.regex + '))', self.sequence)]
        match_pattern= [match.group(1) for match in re.finditer('(?=(' + slim_type_inst.regex + '))', self.sequence)]
        match_results= list(zip(slim_start, match_pattern))
        if len(match_results) >0:
            self.slim_matches_dict[slim_type_inst.slim_id]= []
            for match in match_results:
                if match[0] == 0:
                    pattern= self.sequence[match[0]:match[0] + len(match[1]) + 1]
                    modified_pattern= '-' + pattern[:-1] + str.lower(pattern[-1])
                elif match[0] + len(match[1]) == len(self.sequence):
                    pattern= self.sequence[match[0] - 1:match[0] + len(match[1])]
                    modified_pattern= str.lower(pattern[0]) + pattern[1:]
                else:
                    pattern= self.sequence[match[0] - 1:match[0] + len(match[1]) + 1]
                    modified_pattern= str.lower(pattern[0]) + pattern[1:-1] + str.lower(pattern[-1])
                slim_match_inst= SLiMMatch(dmi_type_inst, slim_type_inst, self, match[0] + 1, match[0] + len(match[1]), modified_pattern)
                # pattern is saved with flanking residues while start and end is saved as realy aa number, i.e. starting from 1
                self.slim_matches_dict[slim_type_inst.slim_id].append(slim_match_inst)

    def read_in_features(self, features_path):
        with open(features_path + '/IUPred_long/' + self.protein_id + '_iupredlong.txt', 'r') as f:
            lines= [line.strip() for line in f.readlines()]
            for line in lines[1:]:
                self.IUPredLong_scores.append(float(line.split('\t')[2]))
        with open(features_path + '/IUPred_short/' + self.protein_id + '_iupredshort.txt', 'r') as f:
            lines= [line.strip() for line in f.readlines()]
            for line in lines[1:]:
                self.IUPredShort_scores.append(float(line.split('\t')[2]))
        with open(features_path + '/Anchor/' + self.protein_id + '_anchor.txt', 'r') as f:
            lines= [line.strip() for line in f.readlines()]
            for line in lines[1:]:
                self.Anchor_scores.append(float(line.split('\t')[2]))
        with open(features_path + '/Domain_overlap/' + self.protein_id + '_domain_overlap.txt', 'r') as f:
            lines= [line.strip() for line in f.readlines()]
            for line in lines[1:]:
                self.DomainOverlap_scores.append(float(line.split('\t')[2]))
        try:
            with open(features_path + '/conservation_scores/' + self.protein_id + '_con.json', 'r') as f:
                data= json.load(f)
            for result in data['Conservation']:
                if 'qfo' in result:
                    self.qfo_RLC_scores= result["qfo"]
                elif 'vertebrates' in result:
                    self.vertebrates_RLC_scores= result["vertebrates"]
                elif 'mammalia' in result:
                    self.mammalia_RLC_scores= result["mammalia"]
                elif 'metazoa' in result:
                    self.metazoa_RLC_scores= result["metazoa"]
        except:
            print(f'{self.protein_id} does not have a conservation score file.')
        if any(self.IUPredLong_scores):
            print(f'{self.protein_id} IUPred long scores saved.')
        if any(self.qfo_RLC_scores):
            print(f'{self.protein_id} qfo RLC scores saved.')
        if any(self.vertebrates_RLC_scores):
            print(f'{self.protein_id} vertebrates RLC scores saved.')
        if any(self.mammalia_RLC_scores):
            print(f'{self.protein_id} mammalia RLC scores saved.')
        if any(self.metazoa_RLC_scores):
            print(f'{self.protein_id} metazoa RLC scores saved.')

    def calculate_features_scores(self):
        for slim_id, slim_match in self.slim_matches_dict.items():
            for slim_match_inst in slim_match:
                slim_match_inst.get_slim_match_features()

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
    def __init__(self, dmi_type_inst, slim_type_inst, prot_inst, start, end, pattern): # start and end numbered starting from 1
        self.dmi_type_inst= dmi_type_inst # saved as DMIType obj because I need domain match info for domain enrichment calculation
        self.slim_type_inst= slim_type_inst # saved as SLiMType object because I need regex from slim_type_inst for features calculation
        self.prot_inst= prot_inst # saved as a protein obj because I need feature attributes from protein for features calculation
        self.start= start
        self.end= end
        self.pattern= pattern
        self.IUPredLong= None
        self.IUPredShort= None
        self.Anchor= None
        self.DomainOverlap= None
        self.qfo_RLC= None
        self.qfo_RLCvar= None
        self.vertebrates_RLC= None
        self.vertebrates_RLCvar= None
        self.mammalia_RLC= None
        self.mammalia_RLCvar= None
        self.metazoa_RLC= None
        self.metazoa_RLCvar= None
        self.DomainEnrichment_pvalue= None
        self.DomainEnrichment_zscore= None
        self.vertex_with_domain_in_real_network= None

    def get_slim_match_features(self):
        defined_positions_url= 'http://slim.icr.ac.uk/restapi/functions/defined_positions?'
        start= int(self.start)
        end= int(self.end)
        pattern= self.pattern
        regex= self.slim_type_inst.regex
        self.IUPredLong= float(sum(self.prot_inst.IUPredLong_scores[start-1:end])/(end - start + 1))
        self.IUPredShort= float(sum(self.prot_inst.IUPredShort_scores[start-1:end])/(end - start + 1))
        self.Anchor= float(sum(self.prot_inst.Anchor_scores[start-1:end])/(end - start + 1))
        self.DomainOverlap= float(sum(self.prot_inst.DomainOverlap_scores[start-1:end])/(end - start + 1))
        print(f'Average IUPred & DomainOverlap scores of {self.slim_type_inst.slim_id} at {self.start}-{self.end} calculated for {self.prot_inst.protein_id}.')
        payload= {'motif':regex, 'sequence':pattern}
        try:
            timeout= 61
            r= requests.get(defined_positions_url, params= payload, timeout= timeout)
            if r.status_code== requests.codes.ok:
                response= r.json()
                defined_positions= [start + (ind - 1) for ind in response['indexes']]
                for i, cons_type in enumerate([self.prot_inst.qfo_RLC_scores, self.prot_inst.vertebrates_RLC_scores, self.prot_inst.mammalia_RLC_scores, self.prot_inst.metazoa_RLC_scores]):
                    if any(cons_type):
                        defined_positions_cons_scores= []
                        for pos, score in cons_type.items():
                            if int(pos) in defined_positions:
                                defined_positions_cons_scores.append(score)
                        if any(defined_positions_cons_scores):
                            pmotif= np.product(defined_positions_cons_scores)
                            lnpmotif= -np.log(pmotif)
                            sigmotif= sc.gammaincc(len(defined_positions_cons_scores), lnpmotif)
                            meanRLCprob= np.mean(defined_positions_cons_scores)
                            varRLCprob= sum([abs(x-meanRLCprob) for x in defined_positions_cons_scores])/len(defined_positions_cons_scores)
                            if i == 0:
                                self.qfo_RLC= sigmotif
                                self.qfo_RLCvar= varRLCprob
                            elif i == 1:
                                self.vertebrates_RLC= sigmotif
                                self.vertebrates_RLCvar= varRLCprob
                            elif i == 2:
                                self.mammalia_RLC= sigmotif
                                self.mammalia_RLCvar= varRLCprob
                            elif i == 3:
                                self.metazoa_RLC= sigmotif
                                self.metazoa_RLCvar= varRLCprob
                        else:
                            self.qfo_RLC= 'Problem with regex or motif'
                            self.qfo_RLCvar= 'Problem with regex or motif'
                            self.vertebrates_RLC= 'Problem with regex or motif'
                            self.vertebrates_RLCvar= 'Problem with regex or motif'
                            self.mammalia_RLC= 'Problem with regex or motif'
                            self.mammalia_RLCvar= 'Problem with regex or motif'
                            self.metazoa_RLC= 'Problem with regex or motif'
                            self.metazoa_RLCvar= 'Problem with regex or motif'
            else:
                print('Bad response code: Check regex and matched pattern')
                self.qfo_RLC= 'Problem with regex or motif'
                self.qfo_RLCvar= 'Problem with regex or motif'
                self.vertebrates_RLC= 'Problem with regex or motif'
                self.vertebrates_RLCvar= 'Problem with regex or motif'
                self.mammalia_RLC= 'Problem with regex or motif'
                self.mammalia_RLCvar= 'Problem with regex or motif'
                self.metazoa_RLC= 'Problem with regex or motif'
                self.metazoa_RLCvar= 'Problem with regex or motif'
        except:
            print('SLiM server not responding')
            self.qfo_RLC= f'SLiM server not responding after {timeout} time-out.'
            self.qfo_RLCvar= f'SLiM server not responding after {timeout} time-out.'
            self.vertebrates_RLC= f'SLiM server not responding after {timeout} time-out.'
            self.vertebrates_RLCvar= f'SLiM server not responding after {timeout} time-out.'
            self.mammalia_RLC= f'SLiM server not responding after {timeout} time-out.'
            self.mammalia_RLCvar= f'SLiM server not responding after {timeout} time-out.'
            self.metazoa_RLC= f'SLiM server not responding after {timeout} time-out.'
            self.metazoa_RLCvar= f'SLiM server not responding after {timeout} time-out.'
        if any(self.prot_inst.networks):
            num_rand_networks= len(self.prot_inst.networks) - 1 # The first network is the real network
            domain_match_list= self.dmi_match_inst.domain_interface_match # returns a list of domain match that is found in DMI match
            vertices_with_overlapping_domains= {} # saved as {network_id:protein_count_with_overlapping_domain}
            # domain_match= set([d.domain_id for d in domain_match_list])
            # for network_id, network in self.prot_inst.networks.items():
            #     count= 0
            #     for partner in network:
            #         partner_domains= set(partner.domain_matches_dict.keys())
            #         if domain_match.intersection(partner_domains) == domain_match: # strict filtering because the partner MUST have the same domain interface i.e. 1 domain or 2 domains in one protein
            #             count += 1
            #     vertices_with_overlapping_domains[int(network_id)]= count
            if len(cognate_domains) == 1: # Cases of 1 domain or 2 domains in one protein
                domains= set(self.dmi_type_inst.domain_interfaces[0].domain_dict.keys())
                for network_id, network in self.prot_inst.networks.items():
                    count= 0
                    for partner in network:
                        partner_domains= set(partner.domain_matches_dict.keys())
                        if domains.intersection(partner_domains) == domains: # strict filtering because the partner MUST have the same domain interface i.e. 1 domain or 2 domains in one protein
                            count += 1
                    vertices_with_overlapping_domains[int(network_id)]= count
            else: # cases of len(cognate_domains) == 2 -> 2 domains in 2 proteins
            # Need to find which domain the interaction partner has, and find only those proteins in the real network
                domains= set()
                for domain_intf_obj in cognate_domains:
                    for domain_id in domain_intf_obj.domain_dict.keys():
                        domains.add(domain_id)
                for network_id, network in self.prot_inst.networks.items():
                    count= 0
                    for partner in network:
                        partner_domains= set(partner.domain_matches_dict.keys())
                        if any(domains.intersection(partner_domains)): # loose filtering because the partner only needs to have at least oen of the domain interface
                            count += 1
                    vertices_with_overlapping_domains[int(network_id)]= count
            count_in_real_network= vertices_with_overlapping_domains[0] # no. of vertices with domain in real network
            self.vertex_with_domain_in_real_network= count_in_real_network
            num_network_more_equal_real= len([v for v in list(vertices_with_overlapping_domains.values())[1:] if v >= count_in_real_network])
            rand_network_mean= np.mean(list(vertices_with_overlapping_domains.values())[1:])
            rand_network_std= np.std(list(vertices_with_overlapping_domains.values())[1:])
            self.DomainEnrichment_pvalue= num_network_more_equal_real/num_rand_networks
            if rand_network_std == 0:
                self.DomainEnrichment_zscore= dummy_value
            else:
                self.DomainEnrichment_zscore= (count_in_real_network - rand_network_mean) / rand_network_std

class DomainType(protein_interaction_interfaces.DomainType):
    def __init__(self, domain_id):
        super().__init__(domain_id)
        self.dmi_types= []

class DomainInterfaceMatch(protein_interaction_interfaces.DomainInterfaceMatch):
    def __init__(self, slim_id, domain_interface, domain_matches):
        super().__init__(domain_interface, domain_matches)
        self.slim_id= slim_id

class DMIType:
    def __init__(self, slim_id):
        self.slim_id= slim_id
        self.domain_interfaces= [] # [{'domain_id1': 1, 'domain_id2': 10}] / [{'domain_id3':1}, {'domain_id4': 4}], [domain_dict1, domain_dict2]

class DMIMatch:
    def __init__(self,slim_protein,domain_protein,slim_match,domain_interface_match):
        self.slim_protein = slim_protein
        self.domain_protein = domain_protein
        self.slim_match = slim_match
        self.domain_interface_match = domain_interface_match
        self.score = None

### Rework this part
class InterfaceHandling(protein_interaction_interfaces.InterfaceHandling):
    def __init__(self, prot_path, slim_type_file, dmi_type_file, smart_domain_types_file, pfam_domain_types_file, smart_domain_matches_json_file, pfam_domain_matches_json_file, features_path, PPI_file=  None, network_path= None):
        super().__init__(prot_path, PPI_file)
        self.slim_types_dict= {}
        self.dmi_types_dict= {}
        self.slim_type_file= slim_type_file
        self.dmi_type_file= dmi_type_file
        self.smart_domain_types_file= smart_domain_types_file
        self.pfam_domain_types_file= pfam_domain_types_file
        self.smart_domain_matches_json_file= smart_domain_matches_json_file
        self.pfam_domain_matches_json_file= pfam_domain_matches_json_file
        self.features_path= features_path
        if network_path != None:
            self.network_path= network_path # for protein network's info when PRS and RRS needs PRS PPI filtered from real network to prevent circularity
        else:
            self.network_path= self.features_path

    def read_in_proteins(self):
        file_names= [file_name for file_name in glob.glob(self.prot_path + '/*')]
        for file_name in file_names:
            with open(file_name,'r') as file:
                lines = [line.strip() for line in file.readlines()]
            for line in lines:
                if line[0] == '>':
                    protein_id= line[1:]
                    prot_inst= Protein(protein_id)
                    self.proteins_dict[protein_id]= prot_inst
                else:
                    self.proteins_dict[protein_id].sequence = line
        print(f'{len(self.proteins_dict)} proteins read in.')

    # def read_in_known_PPIs(self):
    #     with open(self.PPI_file, 'r') as file:
    #         lines= [line.strip() for line in file.readlines()] # PRS saved as .tsv
    #     for line in lines:
    #         tab= line.split('\t')
    #         PPI_instance= sorted(list([tab[0], tab[1]])) # a sorted protein pair as list
    #         self.known_PPIs.append(tuple(PPI_instance)) # PPI pair saved as tuple
    #     print(f'{len(self.known_PPIs)} PPIs read in.')

    def read_in_slim_types(self):

        """ For some reason, the regex that starts with [] python automatically adds a "" around the string, so one needs to modify these manually """

        with open(self.slim_type_file, 'r') as file:
            lines = [line.strip() for line in file.readlines()]
        for line in lines[5:]:
            tab= line.split('\t')
            slim_id= tab[0]
            slim_type_inst= SLiMType(slim_id)
            slim_type_inst.name= tab[1]
            slim_type_inst.regex= tab[4] #.replace('"', '')
            slim_type_inst.probability= str(tab[5])
            self.slim_types_dict[slim_id] = slim_type_inst
        print(f'{len(self.slim_types_dict)} SLiM types read in.')

    def read_in_DMI_types(self):
        with open(self.dmi_type_file, 'r') as file:
            lines= [line.strip() for line in file.readlines()]
        for line in lines[1:]:
            tab= line.split('\t')
            if tab[6] == '1': # Default use 1 == 1
                slim_id= tab[0]
                domain_id= tab[2]
                if domain_id not in self.domain_types_dict:
                    self.domain_types_dict[domain_id] = DomainType(domain_id)
                domain_id2= ''
                domain_count= int(tab[5])
                domain_count2= ''
                DMI_type_inst= DMIType(slim_id)
                if slim_id not in self.dmi_types_dict.keys():
                    self.dmi_types_dict[slim_id] = DMI_type_inst
                if (len(tab) < 11) or (len(tab) > 10 and tab[11] == '') :
                    domain_interface_inst= DomainInterface()
                    domain_interface_inst.domain_dict[domain_id]= domain_count
                    self.dmi_types_dict[slim_id].domain_interfaces.append(domain_interface_inst)
                    self.domain_types_dict[domain_id].dmi_types.append(self.dmi_types_dict[slim_id])
                elif (len(tab) > 7) & (tab[11] == '1'):
                    domain_id2= tab[7]
                    domain_count2= tab[10]
                    domain_interface_inst= DomainInterface()
                    domain_interface_inst.domain_dict[domain_id] = domain_count
                    domain_interface_inst.domain_dict[domain_id2] = domain_count2
                    if domain_id2 not in self.domain_types_dict:
                        self.domain_types_dict[domain_id2] = DomainType(domain_id2)
                    self.dmi_types_dict[slim_id].domain_interfaces.append(domain_interface_inst)
                    self.domain_types_dict[domain_id2].dmi_types.append(self.dmi_types_dict[slim_id])
                elif (len(tab) > 7) & (tab[11] == '0'):
                    continue
        print(f'{len(self.dmi_types_dict)} DMI types read in.')

    def read_in_domain_types(self): # this one reads in only domains involved in DMI
        for domain_types_file in list([self.smart_domain_types_file, self.pfam_domain_types_file]):
            with open(domain_types_file, 'r') as file:
                lines= [line.strip() for line in file.readlines()]
            for line in lines[2:]:
                tab= line.split('\t')
                domain_id= tab[1]
                if domain_id in self.domain_types_dict:
                    self.domain_types_dict[domain_id].name= tab[0]
                    self.domain_types_dict[domain_id].descr= tab[2]
                    self.domain_types_dict[domain_id].DomainFreqbyProtein= tab[3]
                    self.domain_types_dict[domain_id].DomainFreqinProteome= tab[4]
                    if domain_id[:2] == 'PF':
                        self.domain_types_dict[domain_id].source= 'PFAM'
                    elif domain_id[:2] == 'SM':
                        self.domain_types_dict[domain_id].source= 'SMART'
                    print(f'{domain_id} interacts {len(self.domain_types_dict[domain_id].dmi_types)} DMI types.')
            print(f'{domain_types_file} read in.')
        print(f'{len(self.domain_types_dict)} read in.')

    def read_in_domain_matches(self): # this one reads in only domain matches involved in DMI
        if len(self.domain_types_dict) == 0:
            print('Warning: No DMI types have been read in before')
            print('Reading in domain types now...')
            self.read_in_domain_types()
        for domain_matches_json_file in [self.smart_domain_matches_json_file, self.pfam_domain_matches_json_file]:
            with open(domain_matches_json_file, 'r') as f:
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
                    print(f'{protein_id} has {len(self.proteins_dict[protein_id].domain_matches_dict)} domain types match.')

    def read_in_networks(self):
        if self.network_path == self.features_path:
            file_names= [file_name for file_name in glob.glob(self.features_path + '/Protein_networks/*')]
        else:
            file_names= [file_name for file_name in glob.glob(self.network_path + '/*')]
        for file_name in file_names:
            prot_file= file_name.split('/')[-1]
            prot_id= prot_file.split('_')[0]
            with open(file_name,'r') as file:
                lines = [line.strip() for line in file.readlines()]
            for line in lines[1:]:
                tabs= line.split('\t')
                self.proteins_dict[prot_id].networks[int(tabs[0])]= []
                if '|' in tabs[1]:
                    partners= tabs[1].split('|')
                    for partner in partners:
                        self.proteins_dict[prot_id].networks[int(tabs[0])].append(self.proteins_dict[partner]) # network partner saved as prot inst
                else:
                    self.proteins_dict[prot_id].networks[int(tabs[0])].append(self.proteins_dict[tabs[1]])
                self.proteins_dict[prot_id].network_degree= len(self.proteins_dict[prot_id].networks[0])
            print(f'Network file of {prot_id} read in.')

    def create_slim_matches_all_proteins(self):
        for prot_id, prot_inst in self.proteins_dict.items():
            for slim_id, dmi_type_inst in self.dmi_types_dict.items():
                slim_type_inst= self.slim_types_dict[slim_id]
                prot_inst.create_slim_matches(dmi_type_inst, slim_type_inst)
            print(f'{prot_id} has {len(prot_inst.slim_matches_dict.keys())} SLiM types.')

# read_in_iupred_score and calculate_average_iupred can go into Protein class too that takes one protein and run the function using the protein
# and in InterfaceHandling, make another function that runs the Protein.functions globally.
    def read_in_features_all_proteins(self):
        for prot_id, prot_inst in self.proteins_dict.items():
            prot_inst.read_in_features(self.features_path)
        print('SLiM features read in for all proteins.')

    def calculate_features_scores_all_proteins(self): # Dont pre-calculate this for now, only run this function when we want to make a prediction
        for prot_id, prot_inst in self.proteins_dict.items():
            prot_inst.calculate_features_scores()
            print('All SLiM features calculated for all proteins.')

    def find_DMI_matches(self):
        for protpair, protpair_inst in self.protein_pairs_dict.items():
            protein_pairs= [(protpair[0], protpair[1]), (protpair[1], protpair[0])]
            for protein_pair in protein_pairs:
                prot1_domain_matches_dict= self.proteins_dict[protein_pair[0]].domain_matches_dict
                prot2_slim_matches_dict= self.proteins_dict[protein_pair[1]].slim_matches_dict
                unique_slim_ids = set()
                for domain_id in prot1_domain_matches_dict.keys():
                    slim_id_list= [dmi_type.slim_id for dmi_type in self.domain_types_dict[domain_id].dmi_types]
                    slim_id_match_list= set(list(filter(lambda slim_id: slim_id in prot2_slim_matches_dict, slim_id_list))) # Check if protB has the slim_id of a domain in protA
                    unique_slim_ids = unique_slim_ids.union(slim_id_match_list)
                for slim_id_match in unique_slim_ids:
                    domain_interface_list= self.dmi_types_dict[slim_id_match].domain_interfaces
                    for domain_interface in domain_interface_list:
                        cognate_domains= list(domain_interface.domain_dict.keys())
                        if set(cognate_domains).intersection(set(prot1_domain_matches_dict.keys())) == set(cognate_domains):
                            if len(cognate_domains) == 1:
                                domain_interface_match_inst= DomainInterfaceMatch(slim_id_match, domain_interface, [prot1_domain_matches_dict[cognate_domains[0]]])
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
                                print(f'{protein_pair[0]} has {cognate_domains[0]}, {cognate_domains[1]} and {protein_pair[1]} has {slim_id_match}')
                                for slim_match in prot2_slim_matches_dict[slim_id_match]:
                                    DMIMatch_inst= DMIMatch(protein_pair[1], protein_pair[0], slim_match, domain_interface_match_inst)
                                    if slim_id_match not in protpair_inst.dmi_matches_dict.keys():
                                        protpair_inst.dmi_matches_dict[slim_id_match] = []
                                    protpair_inst.dmi_matches_dict[slim_id_match].append(DMIMatch_inst)
        print('DMI matching completed for all proteins.')
