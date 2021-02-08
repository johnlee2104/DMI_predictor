import sys
sys.path.append('/Users/johnlee/Coding/Python/DMI/')
import re, random, json, itertools
from iupred2a import iupred2a_lib
# can be inherited
class SLiMType:

    def __init__(self,slim_id):

        self.slim_id = slim_id # ELM accession number
        self.name = ''
        self.regex = None
        self.taxonomy = None
        self.probability = None


# can be inherited
class DomainType:

    def __init__(self,domain_id):

        self.domain_id = domain_id # = domain hmm
        self.name = '' # domain short name
        self.descr = '' # domain long name
        self.source = None
        self.frequency = None # can get this from UniPort refseq
        # list that contains all the DMI types with that domain type
        self.dmi_types = []


# can be inherited
class DMIType:

    def __init__(self,slim_id):

        self.slim_id = slim_id
        # list of DomainInterface instances -> different instances represent different domain interfaces that are possible for this DMI,
        # one DomainInterface instance represents the constraints on the domain side that form one interface
        self.domain_interfaces = []

    # def __repr__(self):
    #     return '{}'.format(self.domain_interfaces)

# can be inherited
class DomainInterface:

    def __init__(self):

        # {domain_id1:domain_count1,domain_id2:domain_count2}
        self.domain_dict = {}

    # def __repr__(self):
    #     return '{}'.format(self.domain_dict)


# can be inherited
class DMIMatch:

    def __init__(self,slim_protein,domain_protein,slim_match,domain_interface_match):

        self.slim_protein = slim_protein
        self.domain_protein = domain_protein
        self.slim_match = slim_match
        self.domain_interface_match = domain_interface_match
        self.score = None
# KL mentioned that we should add DMIMatch.slim_id but I think that we can get that form slim_match object...

# can be inherited
class SLiMMatch:

    def __init__(self,slim_id,start,end):

        self.slim_id = slim_id
        self.start = start
        self.end = end
#        self.Cons = None
        self.IUPLong = None
        self.IUPShort = None
        self.Anchor = None


# can be inherited
class DomainMatch:

    def __init__(self,domain_id,start,end):

        self.domain_id = domain_id
        self.start = start
        self.end = end
        self.evalue = None


# can be inherited
class DomainInterfaceMatch:

    def __init__(self,slim_id, domain_interface, domain_matches):

        self.slim_id= slim_id
        self.domain_interface = domain_interface
        # should only contain more than one domain match for domain interfaces that either
        # have a domain count > 1 or that consist of more than 1 different domain types
        self.domain_matches = domain_matches # saves the DomainMatch obj JL: saved as tuple if its two domains to form an interface
        self.domain_matches_count = [] # saves an integer that counts the number of a domain match in the protein (Added by John)
#        self.domain_interface_coverage = None
#        self.min_dom_freq = None


# can be inherited
class Protein:

    def __init__(self,protein_id):

        self.protein_id = protein_id
        self.name = ''
        self.descr = ''
        self.sequence = ''
        # {domain_id:[domain_match1,domain_match2,...]}
        self.domain_matches_dict = {}
        # {slim_id1:[slim_match1,slim_match2,...],slim_id2:[slim_match1,slim_match2,...]}
        self.slim_matches_dict = {}


# can be inherited
class ProteinPair:

    def __init__(self,protA,protB):

        self.proteinA = protA
        self.proteinB = protB
        # {slim_id:[dmi_match1,dmi_match2,...]}
        self.dmi_matches_dict = {}


# can be inherited
class RRSFormation:

    def __init__(self,RRS_type):

        # a name for the kind of RRS that we want to create
        self.RRS_type = RRS_type
        self.proteins_dict = {}
        self.SLiM_types_dict = {}
        self.domain_types_dict = {}
        self.DMI_types_dict = {}
        self.known_PPIs = []
        self.protein_pairs_dict = {}
        self.RRS_instances = []

    # function that creates a protein instance for every protein and fills it with ID and seq info
    # fills the dict self.proteins_dict with protein instances: {protein_id:protein_instance}
    # can be inherited
    def read_in_proteins(self,protein_id_seq_file): # will probably store protein_id and protein_seq together

        # file1 = open(protein_id_seq_file,'r')
        # entries = file1.readlines()
        # file1.close()
        # for line in entries:
        #     prot_id = line[:-1]
        #     prot_inst = Protein(prot_id)
        #     self.proteins_dict[prot_id] = prot_inst

        with open(protein_id_seq_file,'r') as file:
            lines = [line.strip() for line in file.readlines()]
        # for each uniprot sequence:
            # get the protein_id -> prot_id
            # get the sequence -> seq
        for line in lines:
            if line[0] == '>':
                tab= line.split('|')
                protein_id= tab[1]
                prot_inst= Protein(protein_id)
                self.proteins_dict[protein_id]= prot_inst
            else:
                self.proteins_dict[protein_id].sequence += line


    # function that creates SLiM_type instances and saves them in self.SLiM_types_dict: {slim_id: slim_instance}
    # can be inherited
    def read_in_slim_types(self,slim_type_file): # suppose that the file is saved as tsv

        """ For some reason, the regex that starts with [] python automatically adds a "" around the string, so one needs to modify these manually """

        with open(slim_type_file, 'r') as file:
            lines = [line.strip() for line in file.readlines()]
        for line in lines[1:]:
            tab= line.split('\t')
            slim_type_inst= SLiMType(tab[0]) # ELM ID
            slim_type_inst.name= tab[1] # SLiM type name
            slim_type_inst.regex= tab[4].replace('"', '') # For some reason, the regex that starts with [] python automatically adds a "" around the string
                                                            #so one needs to modify these before assigning it to the attribute
            slim_type_inst.probability= float(tab[5])
            self.SLiM_types_dict[tab[0]] = slim_type_inst # in this case I assign attribute because putting the instance as a value of the SLiM_types_dict, is it okay?


    # function that creates DMI_type instances and saves them in self.DMI_types_dict
    # can be inherited
    def read_in_DMI_types(self,DMI_file): # saved as tsv

        with open(DMI_file, 'r') as file:
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
                    domain_interface_inst.domain_dict[domain_id]= domain_count # when it encounters same slim with diff domain, it will create a new entry in the dict
                    # self.DMI_types_dict[slim_id].domain_interfaces.append(domain_interface_inst)
                # elif (len(tab) > 7) & (tab[11] == ''):
                #     domain_interface_inst= DomainInterface()
                #     domain_interface_inst.domain_dict[domain_id]= domain_count # when it encounters same slim with diff domain, it will create a new entry in the dict
                #     self.DMI_types_dict[slim_id].domain_interfaces.append(domain_interface_inst)
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

                # else:
                #     self.domain_types_dict[domain_id].dmi_types.append(self.DMI_types_dict[slim_id])
                # else:
                #     domain_id= tab[2]
                #     domain_count= tab[5]
                #     if len(tab) == 7 :
                #         domain_interface_inst= DomainInterface()
                #         domain_interface_inst.domain_dict[domain_id]= domain_count # when it encounters same slim with diff domain, it will create a new entry in the dict
                #         self.DMI_types_dict[slim_id].domain_interfaces.append(domain_interface_inst)
                #     elif (len(tab) > 7) & (tab[11] == ''):
                #         domain_interface_inst= DomainInterface()
                #         domain_interface_inst.domain_dict[domain_id]= domain_count # when it encounters same slim with diff domain, it will create a new entry in the dict
                #         self.DMI_types_dict[slim_id].domain_interfaces.append(domain_interface_inst)
                #     elif (len(tab) > 7) & (tab[11] == '1'):
                #         domain_id2= tab[7]
                #         domain_count2= tab[10]
                #         domain_interface_inst= DomainInterface()
                #         domain_interface_inst.domain_dict[domain_id] = domain_count
                #         domain_interface_inst.domain_dict[domain_id2] = domain_count2
                #         self.DMI_types_dict[slim_id].domain_interfaces.append(domain_interface_inst)


                # domain_id= tab[2]
                # domain_count= tab[5]
                # if len(tab) == 7 :
                #     domain_interface_inst= DomainInterface()
                #     domain_interface_inst.domain_dict[domain_id]= domain_count # when it encounters same slim with diff domain, it will create a new entry in the dict
                #     self.DMI_types_dict[slim_id].domain_interfaces.append(domain_interface_inst)
                # elif (len(tab) > 7) & (tab[11] == ''):
                #     domain_interface_inst= DomainInterface()
                #     domain_interface_inst.domain_dict[domain_id]= domain_count # when it encounters same slim with diff domain, it will create a new entry in the dict
                #     self.DMI_types_dict[slim_id].domain_interfaces.append(domain_interface_inst)
                # elif (len(tab) > 7) & (tab[11] != ''):
                #     domain_id2= tab[7]
                #     domain_count2= tab[10]
                #     domain_interface_inst= DomainInterface()
                #     domain_interface_inst.domain_dict[domain_id] = domain_count
                #     domain_interface_inst.domain_dict[domain_id2] = domain_count2
                #     self.DMI_types_dict[slim_id].domain_interfaces.append(domain_interface_inst)

    # function that creates domain_type instances for those domains that occur in the DMI types
    # it saves the domain_type instances in the self.domain_types_dict: {domain_id: domain_type instance}
    # can be inherited
    # def read_in_domain_types(self,domain_types_file): # Using the DMI list I annotated
    #
    #     with open(domain_types_file, 'r') as file: # not sure why make the distinction between Pfam and SMART, I can just read the whole list in and use the .source to specify source
    #         lines= [line.strip() for line in file.readlines()] #KL: save only domains with default use == 1, remember to save frequency attr
    #     for line in lines[1:]:
    #         tab= line.split('\t')
    #         if len(tab) == 7: # this is to deal with the fact that some rows have comments at the very end and these rows will have many empty tabs in between
    #             if tab[6] == '1': # simple cases of just one cognate domain
    #                 slim_id= tab[0]
    #                 domain_id= tab[2] # domain HMM
    #                 if domain_id not in self.domain_types_dict.keys():
    #                     domain_type_inst = DomainType(domain_id)
    #                     domain_type_inst.name= tab[3]
    #                     domain_type_inst.descr= tab[4]
    #                     if domain_id[:2] == 'PF':
    #                         domain_type_inst.source= 'PFAM'
    #                     else:
    #                         domain_type_inst.source= 'SMART' # pfam/SMART
    #                     domain_type_inst.dmi_types.append(self.DMI_types_dict[slim_id])
    #                     self.domain_types_dict[domain_id] = domain_type_inst
    #                 elif domain_id in self.domain_types_dict.keys():
    #                     self.domain_types_dict[domain_id].dmi_types.append(self.DMI_types_dict[slim_id])
    #         else: # these rows have either one cognate domain and one comment at the very end, or two domains in the same row
    #             if (tab[6] == '1') & (tab[11] == ''): # simple cases with one cognate domain and a comment at the very end
    #                 slim_id= tab[0]
    #                 domain_id= tab[2] # domain HMM
    #                 if domain_id not in self.domain_types_dict.keys():
    #                     domain_type_inst = DomainType(domain_id)
    #                     domain_type_inst.name= tab[3]
    #                     domain_type_inst.descr= tab[4]
    #                     if domain_id[:2] == 'PF':
    #                         domain_type_inst.source= 'PFAM'
    #                     else:
    #                         domain_type_inst.source= 'SMART' # pfam/SMART
    #                     domain_type_inst.dmi_types.append(self.DMI_types_dict[slim_id])
    #                     self.domain_types_dict[domain_id] = domain_type_inst
    #                 elif domain_id in self.domain_types_dict.keys():
    #                     self.domain_types_dict[domain_id].dmi_types.append(self.DMI_types_dict[slim_id])
    #             if (tab[6] == '1') & (tab[11] == '1'): # special cases with two cognate domains in the same row
    #                 slim_id= tab[0]
    #                 domain_id= tab[2] # domain HMM
    #                 if domain_id not in self.domain_types_dict.keys():
    #                     domain_type_inst = DomainType(domain_id)
    #                     domain_type_inst.name= tab[3]
    #                     domain_type_inst.descr= tab[4]
    #                     if domain_id[:2] == 'PF':
    #                         domain_type_inst.source= 'PFAM'
    #                     else:
    #                         domain_type_inst.source= 'SMART' # pfam/SMART
    #                     domain_type_inst.dmi_types.append(self.DMI_types_dict[slim_id])
    #                     self.domain_types_dict[domain_id] = domain_type_inst
    #                 elif domain_id in self.domain_types_dict.keys():
    #                     self.domain_types_dict[domain_id].dmi_types.append(self.DMI_types_dict[slim_id])
    #
    #                 domain_id2= tab[7] # second cognate domain in the same row
    #                 if domain_id2 not in self.domain_types_dict.keys():
    #                     domain_type_inst = DomainType(domain_id2)
    #                     domain_type_inst.name= tab[8]
    #                     domain_type_inst.descr= tab[9]
    #                     if domain_id2[:2] == 'PF':
    #                         domain_type_inst.source= 'PFAM'
    #                     else:
    #                         domain_type_inst.source= 'SMART' # pfam/SMART
    #                     domain_type_inst.dmi_types.append(self.DMI_types_dict[slim_id])
    #                     self.domain_types_dict[domain_id2] = domain_type_inst
    #                 elif domain_id2 in self.domain_types_dict.keys():
    #                     self.domain_types_dict[domain_id2].dmi_types.append(self.DMI_types_dict[slim_id])


    def read_in_domain_types(self, domain_types_file):

        with open(domain_types_file, 'r') as file:
            lines= [line.strip() for line in file.readlines()]
        for line in lines[2:]:
            tab= line.split('\t')
            name= tab[0]
            domain_id= tab[1]
            if len(tab) > 2:
                descr= tab[2]
            if domain_id not in self.domain_types_dict: # restrict the domain match to only those in DMI types list
                continue
            self.domain_types_dict[domain_id].name= name
            self.domain_types_dict[domain_id].descr= descr
            if domain_id[:2] == 'PF':
                self.domain_types_dict[domain_id].source= 'PFAM'
            elif domain_id[:2] == 'SM':
                self.domain_types_dict[domain_id].source= 'SMART'

    # function that generates slim_match instances for every protein in self.proteins_dict and saves them in
    # the protein instances
    # can be inherited
    def create_slim_matches(self):

        for prot_id, prot_inst in self.proteins_dict.items():
            for slim_id, slim_type_inst in self.SLiM_types_dict.items():
                slim_start= [match.start() + 1 for match in re.finditer('(?=(' + slim_type_inst.regex + '))', prot_inst.sequence)]
                match_pattern= [match.group(1) for match in re.finditer('(?=(' + slim_type_inst.regex + '))', prot_inst.sequence)]
                match_results= list(zip(slim_start, match_pattern))
                if len(match_results) >0:
                    print(prot_id + ' has ' + str(len(match_results)) + ' matches of ' + slim_id + ' ' + slim_type_inst.name)
                    prot_inst.slim_matches_dict[slim_id]= []
                    for match in match_results:
                        slim_match_inst= SLiMMatch(slim_id, match[0], match[0] + len(match[1]) - 1)
                        prot_inst.slim_matches_dict[slim_id].append(slim_match_inst)

    def calculate_slim_matches_iupred(self):

        for prot_id, prot_inst in self.proteins_dict.items():
            seq= prot_inst.sequence
            result= iupred2a_lib.iupred(seq)
            result2= iupred2a_lib.iupred(seq, mode= 'short')
            for slim_id, slim_match in prot_inst.slim_matches_dict.items():
                for slim_match_inst in slim_match:
                    start= int(slim_match_inst.start)
                    end= int(slim_match_inst.end)
                    slim_match_inst.IUPLong= float(sum(result[start-1:end])/(end - start + 1))
                    slim_match_inst.IUPShort= float(sum(result2[start-1:end])/(end - start + 1))


    def calculate_slim_matches_anchor(self):

        for prot_id, prot_inst in self.proteins_dict.items():
            seq= prot_inst.sequence
            result= iupred2a_lib.anchor2(seq)
            for slim_id, slim_match in prot_inst.slim_matches_dict.items():
                for slim_match_inst in slim_match:
                    start= int(slim_match_inst.start)
                    end= int(slim_match_inst.end)
                    slim_match_inst.Anchor= float(sum(result[start-1:end])/(end - start + 1))

    # function that reads in the domain matches from the SMART output file for the domain types in the
    # self.domain_types_dict and the proteins that are in self.proteins_dict
    # it creates domain match instances that are saved in the protein instances
    # can be inherited
    def read_in_domain_matches(self,json_file):

        with open(json_file) as f:
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
        #                 if prot not in all_domain_matches_dict:
        #                     all_domain_matches_dict[prot]= []
        #                 all_domain_matches_dict[prot].append((domain, start, end, source, evalue))
        # for prot_id in self.proteins_dict.keys(): # KL: save the domain matches right away, without making new dictionary
        #     if prot_id in all_domain_matches_dict.keys(): # check if the protein_id has a domain match
        #         domain_matches= all_domain_matches_dict[prot_id]
        #         for domain_match in domain_matches:
        #             for domain_id, domain_type_inst in self.domain_types_dict.items(): # Save evalue to the attr
        #                 if (domain_match[0] == domain_type_inst.name) & (domain_match[3] == domain_type_inst.source): # name same and source same
        #                     domain_match_inst= DomainMatch(domain_type_inst, domain_match[1], domain_match[2])
        #                     domain_match_inst.evalue= domain_match[4]
        #                     if domain_id not in self.proteins_dict[prot_id].domain_matches_dict.keys():
        #                         self.proteins_dict[prot_id].domain_matches_dict[domain_id]= []
        #                         self.proteins_dict[prot_id].domain_matches_dict[domain_id].append(domain_match_inst)
        #                     else:
        #                         self.proteins_dict[prot_id].domain_matches_dict[domain_id].append(domain_match_inst)


    # function that reads in all known PPIs by which the random pairs will need to be filtered
    # saves the PPIs in self.known_PPIs list: [(protA,protB),(protC,protD),...] with protA < protB and
    # protC < protD
    # JL: remember to read in two PPI files, one from PRS and one from IntAct
    # can be inherited
    def read_in_known_PPIs(self,PPI_file):

        with open(PPI_file, 'r') as file:
            lines= [line.strip() for line in file.readlines()] # PRS saved as .tsv
        for line in lines:
            tab= line.split('\t')
            PPI_instance= sorted(list([tab[0], tab[1]])) # a sorted protein pair as list
            self.known_PPIs.append(tuple(PPI_instance)) # PPI pair saved as tuple


    # function that creates a given number of random protein pairs from the proteins in self.proteins_dict
    # and saves them in self.protein_pairs_dict: {(protA,protB): protein_pair instance} where protA < protB
    # when sorted by name
    # checks for every random pair whether it is a known PPI in which case it will be rejected
    # specific function for this module
    def make_random_protein_pairs(self):

        pairs_set= set()

        for p in itertools.combinations(self.proteins_dict.keys(), 2):
            random_protein_pair= sorted(p)
            pairs_set.add(tuple(random_protein_pair))
        for protein_pair in pairs_set.difference(set(self.known_PPIs)):
            if protein_pair not in self.protein_pairs_dict:
                self.protein_pairs_dict[protein_pair] = ProteinPair(protein_pair[0], protein_pair[1])

        # number_of_pairs= 0
        #
        # while True:
        #     random_protein_pair= sorted(random.sample(self.proteins_dict.keys(), 2))
        #     if tuple(random_protein_pair) not in self.known_PPIs:
        #         self.protein_pairs_dict[tuple(random_protein_pair)] = ProteinPair(random_protein_pair[0], random_protein_pair[1])
        #         number_of_pairs += 1
        #         if number_of_pairs == num_pairs:
        #             break

    # needs to find all DMI matches for protein pair in either orientation: protA has SLiM and protB has
    # domain and vice versa
    # for each protein pair:
        # take protA as domain protein
        # for each domain type that has at least 1 match in that protein
            # for each DMI type with that domain type:
                # test whether suitable SLiM match in protB
                    # for each SLiM match and domain match of this DMI type -> save dmi match in
                    # protein_pair.dmi_matches_dict
        # take protB as domain protein
        # repeat all the previous steps
    # can be inherited
    def find_DMI_matches(self):

        # for protpair, protpair_inst in self.protein_pairs_dict.items():
        #     prot1_domain_matches_dict= self.proteins_dict[protpair[0]].domain_matches_dict
        #     prot2_slim_matches_dict= self.proteins_dict[protpair[1]].slim_matches_dict
        #     for domain_id in prot1_domain_matches_dict.keys():
        #         dmi_types= self.domain_types_dict[domain_id].dmi_types # saved as a list of dmi types object from the class DMIType
        #         for dmi_type in dmi_types:
        #             slim_type= dmi_type.slim_type #returns slim_type as literal string
        #             domain_interfaces= self.DMI_types_dict[slim_type].domain_interfaces # return as a list of domain_interface object
        #             for slim_id, slim_type_inst in self.SLiM_types_dict.items():
        #                 if slim_type_inst.name == slim_type:
        #                     if slim_id in prot2_slim_matches_dict.keys():
        #                         for slim_match_inst in prot2_slim_matches_dict[slim_id]:
        #                             DMIMatch_inst= DMIMatch(protpair[1], protpair[2], slim_match_inst, )
        #             if dmi_type.slim_type in prot2_slim_matches_dict.keys()
        #             domain_interfaces= self.DMI_types_dict[slim_type].domain_interfaces
        #             cognate_domain_dicts= [domain_interface.domain_dict for domain_interface in domain_interfaces]

        # JL: my idea is that when there's a dmi match in a random pair of protein, I will save the domain match as an object of class DomainInterfaceMatch, under the .domain_interface_match attr of class DMIMatch
        # and I will save the domain match as an object under the class DomainMatch, in the .domain_match attr of DomainInterfaceMatch. I will then count how many of the same domain match is found in the domain_match
        # of the protein and the number of domain match in the .domain_matches_count attr

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
                                domain_interface_match_inst.domain_matches_count.append(len(prot1_domain_matches_dict[cognate_domains[0]]))
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




        #         for slim_id in prot2_slim_matches_dict.keys():
        #             domain_interfaces= self.DMI_types_dict[slim_id].domain_interfaces # returns a list
        #             for domain_interface in domain_interfaces:
        #                 cognate_domain_dict= domain_interface.domain_dict
        #                 if len(cognate_domain_dict) == 1: # simple cases, trimeric complex or dimeric with 1 slim 1 domain
        #                     for domain_match, domain_match_inst in prot1_domain_matches_dict.items():
        #                         if domain_match in cognate_domain_dict.keys():
        #                             domain_interface_match_inst= DomainInterfaceMatch(domain_interface, domain_match_inst)
        #                             domain_interface_match_inst.domain_matches_count= len(prot1_domain_matches_dict[domain_match])
        #                             if domain_interface_match_inst.domain_matches_count > 1:
        #                                 print('slim ' + protpair[1] + ' 1 domain ' + str(domain_interface_match_inst.domain_matches_count) + ' ' + protpair[0], protpair, slim_id)
        #                             for slim_match in self.proteins_dict[protpair[1]].slim_matches_dict[slim_id]:
        #                                 DMIMatch_inst= DMIMatch(protpair[1], protpair[0], slim_match, domain_interface_match_inst)
        #                                 if slim_id not in protpair_inst.dmi_matches_dict.keys():
        #                                     protpair_inst.dmi_matches_dict[slim_id] = []
        #                                     protpair_inst.dmi_matches_dict[slim_id].append(DMIMatch_inst)
        #                                 else:
        #                                     protpair_inst.dmi_matches_dict[slim_id].append(DMIMatch_inst)
        #
        #                 elif len(cognate_domain_dict) == 2: # special cases, dimeric complex with 1 slim 2 domains in the same protein
        #                     cognate_domains= list(cognate_domain_dict.keys())
        #                     if set(cognate_domain_dict.keys()).intersection(set(prot1_domain_matches_dict.keys())) == set(cognate_domain_dict.keys()):
        #                         domain_match1= prot1_domain_matches_dict[cognate_domains[0]] # a list
        #                         domain_match2= prot1_domain_matches_dict[cognate_domains[1]] # a list
        #                         domain_interface_match_inst= DomainInterfaceMatch(domain_interface, tuple([domain_match1, domain_match2]))
        #                         domain_interface_match_inst.domain_matches_count= [(cognate_domains[0], len(prot1_domain_matches_dict[cognate_domains[0]])), (cognate_domains[1], len(prot1_domain_matches_dict[cognate_domains[1]]))]
        #                         print('Special case: ' + protpair[1] + ' slim protein ' + protpair[0] + ' domain protein.', protpair, slim_id)
        #                         print(cognate_domain_dict, domain_interface_match_inst.domain_matches_count, slim_id)
        #                         for slim_match in self.proteins_dict[protpair[1]].slim_matches_dict[slim_id]:
        #                             DMIMatch_inst= DMIMatch(protpair[1], protpair[0], slim_match, domain_interface_match_inst)
        #                             if slim_id not in protpair_inst.dmi_matches_dict.keys():
        #                                 protpair_inst.dmi_matches_dict[slim_id] = []
        #                                 protpair_inst.dmi_matches_dict[slim_id].append(DMIMatch_inst)
        #                             else:
        #                                 protpair_inst.dmi_matches_dict[slim_id].append(DMIMatch_inst)
        #
        # # Same code, just changing prot[0] to be slim protein and prot[1] to be domain protein
        # for protpair, protpair_inst in self.protein_pairs_dict.items():
        #     prot1_domain_matches_dict= self.proteins_dict[protpair[1]].domain_matches_dict
        #     prot2_slim_matches_dict= self.proteins_dict[protpair[0]].slim_matches_dict
        #     for slim_id in prot2_slim_matches_dict.keys():
        #         # DMI_type= self.DMI_types_dict[slim_id]
        #         domain_interfaces= self.DMI_types_dict[slim_id].domain_interfaces # returns a list
        #         for domain_interface in domain_interfaces:
        #             cognate_domain_dict= domain_interface.domain_dict
        #             if len(cognate_domain_dict) == 1: # simple cases, trimeric complex or dimeric with 1 slim 1 domain
        #                 for domain_match, domain_match_inst in prot1_domain_matches_dict.items():
        #                     if domain_match in cognate_domain_dict.keys():
        #                         domain_interface_match_inst= DomainInterfaceMatch(domain_interface, domain_match_inst)
        #                         domain_interface_match_inst.domain_matches_count= len(prot1_domain_matches_dict[domain_match])
        #                         if domain_interface_match_inst.domain_matches_count > 1:
        #                             print('slim ' + protpair[0] + ' 1 domain ' + str(domain_interface_match_inst.domain_matches_count) + ' ' + protpair[1], protpair, slim_id)
        #                         for slim_match in self.proteins_dict[protpair[0]].slim_matches_dict[slim_id]:
        #                             DMIMatch_inst= DMIMatch(protpair[0], protpair[1], slim_match, domain_interface_match_inst)
        #                             if slim_id not in protpair_inst.dmi_matches_dict.keys():
        #                                 protpair_inst.dmi_matches_dict[slim_id] = []
        #                                 protpair_inst.dmi_matches_dict[slim_id].append(DMIMatch_inst)
        #                             else:
        #                                 protpair_inst.dmi_matches_dict[slim_id].append(DMIMatch_inst)
        #
        #             elif len(cognate_domain_dict) == 2: # special cases, dimeric complex with 1 slim 2 domains in the same protein
        #                 cognate_domains= list(cognate_domain_dict.keys())
        #                 if set(cognate_domain_dict.keys()).intersection(set(prot1_domain_matches_dict.keys())) == set(cognate_domain_dict.keys()):
        #                     domain_match1= prot1_domain_matches_dict[cognate_domains[0]] # a list
        #                     domain_match2= prot1_domain_matches_dict[cognate_domains[1]] # a list
        #                     domain_interface_match_inst= DomainInterfaceMatch(domain_interface, tuple([domain_match1, domain_match2]))
        #                     domain_interface_match_inst.domain_matches_count= [(cognate_domains[0], len(prot1_domain_matches_dict[cognate_domains[0]])), (cognate_domains[1], len(prot1_domain_matches_dict[cognate_domains[1]]))]
        #                     print('Special case: ' + protpair[0] + ' slim protein ' + protpair[1] + ' domain protein.', protpair, slim_id)
        #                     print(cognate_domains, domain_interface_match_inst.domain_matches_count, slim_id)
        #                     for slim_match in self.proteins_dict[protpair[0]].slim_matches_dict[slim_id]:
        #                         DMIMatch_inst= DMIMatch(protpair[0], protpair[1], slim_match, domain_interface_match_inst)
        #                         if slim_id not in protpair_inst.dmi_matches_dict.keys():
        #                             protpair_inst.dmi_matches_dict[slim_id] = []
        #                             protpair_inst.dmi_matches_dict[slim_id].append(DMIMatch_inst)
        #                         else:
        #                             protpair_inst.dmi_matches_dict[slim_id].append(DMIMatch_inst)

                        # for domain_match in prot1_domain_matches_dict.keys():
                        #     for cognate_domain_dict in cognate_domain_dicts:
                        #         if domain_match in cognate_domain_dict.keys():
                        #             domain_interface_match_inst= DomainInterfaceMatch(domain_interface_inst, prot1_domain_matches_dict[domain_match])
                        #             for slim_match_inst in self.proteins_dict[protpair[1]].slim_matches_dict[slim_id]:
                        #                 DMIMatch_inst= DMIMatch(protpair[1], protpair[0], slim_match_inst, domain_interface_match_inst)
                        #                 if slim_id not in protpair_inst.dmi_matches_dict.keys():
                        #                     protpair_inst.dmi_matches_dict[slim_id] = []
                        #                     protpair_inst.dmi_matches_dict[slim_id].append(DMIMatch_inst)
                        #                 else:
                        #                     protpair_inst.dmi_matches_dict[slim_id].append(DMIMatch_inst)

        # for protpair, protpair_inst in self.protein_pairs_dict.items():
        #     prot1_domain_matches_dict= self.proteins_dict[protpair[1]].domain_matches_dict
        #     prot2_slim_matches_dict= self.proteins_dict[protpair[0]].slim_matches_dict
        #     for slim_id in prot2_slim_matches_dict.keys():
        #         DMI_type_inst= self.DMI_types_dict[slim_id]
        #         domain_interface_inst= DMI_type_inst.domain_interfaces # returns a list of domain_interface instance(s)
        #         cognate_domain_dicts= [domain_interface.domain_dict for domain_interface in domain_interface_inst] # store the domain_dict of each domain_interface instance
        #         for domain_match in prot1_domain_matches_dict.keys():
        #             for cognate_domain_dict in cognate_domain_dicts:
        #                 if domain_match in cognate_domain_dict.keys():
        #                     domain_interface_match_inst= DomainInterfaceMatch(domain_interface_inst, prot1_domain_matches_dict[domain_match])
        #                     for slim_match_inst in self.proteins_dict[protpair[0]].slim_matches_dict[slim_id]:
        #                         DMIMatch_inst= DMIMatch(protpair[0], protpair[1], slim_match_inst, domain_interface_match_inst)
        #                         if slim_id not in protpair_inst.dmi_matches_dict.keys():
        #                             protpair_inst.dmi_matches_dict[slim_id] = []
        #                             protpair_inst.dmi_matches_dict[slim_id].append(DMIMatch_inst)
        #                         else:
        #                             protpair_inst.dmi_matches_dict[slim_id].append(DMIMatch_inst)


    # function that selects out of all random protein pairs and identified DMI matches randomly some to
    # be used in our RRS according to our defined rules
    # for every DMI type, select X number of DMI matches whereby each DMI match should come from a different
    # protein pair for a given DMI type (for other DMI types, the protein pair can occur again)
    # save the selected DMI matches in the list self.RRS_instances
    # specific function for this module
    def select_RRS_instances(self,num_instances):

        # for each DMI type:
            # find all protein pairs that have at least 1 dmi match of this DMI type -> save these
            # protein pairs in a list
            # randomly select X number of protein pairs from this list (without replacement)
            # for each randomly selected protein pair:
                # randomly select a DMI match out of all DMI matches for that protein pair from that DMI type
                # save the selected DMI match in the list self.RRS_instances
        for slim_id in self.DMI_types_dict.keys():
            dmi_matches= []
            for protpair, protpair_inst in self.protein_pairs_dict.items():
                if slim_id in protpair_inst.dmi_matches_dict.keys():
                    dmi_match_insts= protpair_inst.dmi_matches_dict[slim_id] # save all dmi match instances in a list
                    for random_inst in random.sample(dmi_match_insts, 1): # here I randomly select just 1 dmi match inst from the protein pair
                        dmi_matches.append(random_inst)             # to be added into the dmi_matches list, thus restricting 1 dmi match inst per PP
            if len(dmi_matches) <= num_instances:
                self.RRS_instances = self.RRS_instances + dmi_matches
            else:
                selected_dmi_matches = random.sample(dmi_matches, num_instances) # the function sample itself is done without replacement
                self.RRS_instances = self.RRS_instances + selected_dmi_matches

        # return self.RRS_instances


    # function that can be called at any stage in the program to dump the content of the current
    # RRSFormation instance using the pickle module of python and the pickle.dump() function
    # this might be handy when debugging the code and you don't want to start all over reading
    # everything in, or you test the code in jupyter until it works and then you copy-paste it
    # back over here, I heard, there is also ways to do this with Atom
    def dump_RRSFormation_instance(self,outfile):

        from pickle import dump
        with open(outfile, 'wb') as pickle_file: # 'with' automatically closes the file
            dump(self, pickle_file)

    # function that writes out all selected RRS instances into a tab-separated file
    def write_out_RRS_instances(self,outfile):

        file= open(outfile, 'w')
        file.write('\t'.join((' ', 'Elm', 'interactorElm', 'ElmMatch', 'IUPredLong', 'IUPredShort', 'Anchor',
        'interactorDomain', 'Domain_ID1', 'DomainMatch1', 'Evalue1', 'Domain_ID2', 'DomainMatch2', 'Evalue2')))
        file.write('\n')
        for i , inst in enumerate(self.RRS_instances):
            file.write('\t'.join((str(i+1), self.SLiM_types_dict[inst.slim_match.slim_id].name, inst.slim_protein ,
            '-'.join([str(inst.slim_match.start) , str(inst.slim_match.end)]), str(inst.slim_match.IUPLong),
            str(inst.slim_match.IUPShort), str(inst.slim_match.Anchor), inst.domain_protein)))
            for domain_match_list in inst.domain_interface_match.domain_matches:
                domain_id= domain_match_list[0].domain_id
                start= [domain_match.start for domain_match in domain_match_list]
                end= [domain_match.end for domain_match in domain_match_list]
                evalue= [domain_match.evalue for domain_match in domain_match_list]
                match= '|'.join([f"{i[0]}-{i[1]}" for i in zip(start, end)])
                evalues= '|'.join([str(ev) for ev in evalue])
                file.write('\t')
                file.write('\t'.join((domain_id, match, evalues)))
            file.write('\n')
            # elif len(inst.domain_interface_match.domain_interface.domain_dict) == 2:
            #     print('tuple!')
            #     file.write('\t'.join((str(i+1) , inst.slim_match.slim_id.name, inst.slim_protein , str(inst.slim_match.start) , str(inst.slim_match.end)
            #     , inst.domain_protein , list(inst.domain_interface_match.domain_interface.domain_dict.keys())[0] , str(inst.domain_interface_match.domain_match[0].start)
            #     , str(inst.domain_interface_match.domain_match[0].end) , str(inst.domain_interface_match.domain_matches_count[0][1])
            #     , list(inst.domain_interface_match.domain_interface.domain_dict.keys())[1] , str(inst.domain_interface_match.domain_match[1].start)
            #     , str(inst.domain_interface_match.domain_match[1].end), str(inst.domain_interface_match.domain_matches_count[1][1]))))
            #     file.write('\n')
        file.close()
# if __main__ == '__name__':
#
#     lig_pdz_1 = SLiMType('ELM95757585')
#     lig_pdz_1.name = 'LIG_PDZ_1'
#     lig_pdz_1.regex = '.[TS].[VL]'
#
#     rrs_inst = RRSFormation('RRSv1')
#     rrs_inst.read_in_proteins('/Users/luck/uniprot_id_file.txt','/Users/luck/uniprot_seq_file.fasta')
