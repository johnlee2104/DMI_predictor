# This script includes basal objects and functions useful for the prediction of protein interaction interface. These basal objects and functions form the basis for DMI and DDI predictions.
# Author: Chop Yan Lee

import glob, json

class Protein:
    """
    Represents proteins

    Attributes:
        protein_id (str): UniProt ID of a protein
        sequence (str): Amino acid sequence of a protein
        domain_matches_dict (dict): Dictionary that saves domain IDs (Pfam or SMART ID) as key and lists of their DomainMatch objects as value
    """

    def __init__(self, protein_id):
        """
        Initialize a new instance of Protein class

        Args:
            protein_id (str): UniProt ID of the protein
        """
        self.protein_id= protein_id
        self.sequence= ''
        self.domain_matches_dict= {} # {'domain_id': [domain_match1, domain_match2, ...]}

class ProteinPair:
    """
    Represents protein pairs

    Attributes:
        proteinA (str): UniProt ID of the protein A in protein pair
        proteinB (str): UniProt ID of the protein B in protein pair
    """

    def __init__(self, proteinA, proteinB):
        """
        Initialize a new instance of ProteinPair class

        Args:
            proteinA (str): UniProt ID of the protein A
            proteinB (str): UniProt ID of the protein B
        """
        self.proteinA= proteinA
        self.proteinB= proteinB

class DomainType:
    """
    Represents domain type

    Attributes:
        domain_id (str): Pfam or SMART ID of domain type
        name (str): Name of domain type
        source (str): Database of the domain type, i.e. Pfam or SMART
        DomainFreqbyProtein (float): Frequency of the domain type calculated as the fraction of human SwissProt proteins having at least one match of the HMM of the domain type
        DomainFreqinProteome (float): Frequency of the domain type calculated through the frequency of the HMM of the domain type matching in the human SwissProt proteins
    """

    def __init__(self, domain_id):
        """
        Initialize a new instance of DomainType class

        Args:
            domain_id (str): Pfam or SMART ID of domain type
        """
        self.domain_id= domain_id
        self.name= ''
        self.source= ''
        self.DomainFreqbyProtein= None
        self.DomainFreqinProteome= None

class DomainMatch:
    """
    Represents domain match

    Attributes:
        domain_id (str): Pfam or SMART ID of domain type
        start (int): Start of the domain match in a given protein sequence
        end (int): Start of the domain match in a given protein sequence
        evalue (float): Evalue of the domain match
    """

    def __init__(self, domain_id, start, end):
        """
        Initialize a new instance of DomainMatch class

        Args:
            domain_id (str): Pfam or SMART ID of domain type
            start (int): Start of the domain match in a given protein sequence
            end (int): Start of the domain match in a given protein sequence
        """
        self.domain_id= domain_id
        self.start= start
        self.end= end
        self.evalue= None

class DomainInterface:
    """
    Represents the number of a given domain type matching in a protein sequence. Used to set up DMI class where some DMI requires a certain number of tandem repeat occurring in a protein to bind to a motif. For example, some motifs bind to beta barrel structure formed by 6 repeats of WD40 repeat, while others bind to that formed by 7 repeats of WD40 repeat.

    Attributes:
        domain_dict (dict): Dict of domain id as key and count as value to form a domain interface where a given motif type can bind to
    """
    def __init__(self):
        """
        Instantiate DomainInterface 
        """
        self.domain_dict= {}

class DomainInterfaceMatch:
    """
    Represents the domain interface match of a protein that is used for the pairing with a slim_id to form DMI. This class will be inherited in the DMIDB script to form the DMI-specific domain interface.

    Attributes:
        domain_interface (DomainInterface): An instance of DomainInterface
        domain_matches (list of DomainMatch): A list of DomainMatch instances that corresponds to the domain matches found in a given protein
    """

    def __init__(self, domain_interface, domain_matches):
        """
        Instantiate DomainInterfaceMatch

        Args:
            domain_interface (DomainInterface): An instance of DomainInterface
            domain_matches (list of DomainMatch): A list of DomainMatch instances that corresponds to the domain matches found in a given protein
        """
        self.domain_interface= domain_interface
        self.domain_matches= domain_matches # saved as []

class InterfaceHandling:
    """
    A composite class that handles interface information by reading and storing them in this class using the classes that are set up above. This class additionally manipulates the stored information using its own methods or the methods of the classes.

    Attributes:
        PPI_file (str): Absolute path to a file where known PPIs are stored, default None. Used in different RRS_formation script where randomized protein pairs need to checked with known PPIs in order to not accidentally use an interacting pair of proteins as random protein pair
        proteins_dict (dict): Dict that stores uniprot_id as key and an instance of Protein class as value
        domain_types_dict (dict): Dict that stores domain id as key and an instance of DomainType class as value
        known_PPIs (set of tuples): A set that stores known PPI in tuples of protein pairs
        protein_pairs_dict (dict): Dict that stores the uniprot ids of a protein pair in tuple as keys and their ProteinPair instance as value
        prot_path (str): Absolute path to the directory where txt files containing uniprot ids and their sequences are stored
    """

    def __init__(self, prot_path, PPI_file= None):
        """
        Instantiate InterfaceHandling

        Args:
            prot_path (str): Absolute path to the directory where txt files containing uniprot ids and their sequences are stored
            PPI_file (str): Absolute path to a file where known PPIs are stored, default None. Used in different RRS_formation script where randomized protein pairs need to checked with known PPIs in order to not accidentally use an interacting pair of proteins as random protein pair
        """
        if PPI_file != None:
            self.PPI_file= PPI_file
        self.proteins_dict= {}
        self.domain_types_dict= {}
        self.known_PPIs= set()
        self.protein_pairs_dict= {}
        self.prot_path= prot_path

    def read_in_proteins(self, only_canonical= True):
        """
        Look into self.prot_path and reads in proteins and their sequences using the Protein class and store them in self.proteins_dict

        Args:
            only_canonical (bool): True if only use canonical sequences, else all sequences
        """
        if only_canonical == False:
            file_names= [file_name for file_name in glob.glob(self.prot_path + '/*')]
        else: # only_canonical == True
            file_names= [file_name for file_name in glob.glob(self.prot_path + '/*') if '-' not in file_name]
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

    def read_in_domain_types(self, domain_types_file): # This one reads in all domain types, useful for DDI predictor
        """
        Reads domain_types_file that contains domain type information and sets up DomainType instances that are subsequently stored in self.domain_types_dict

        Args:
            domain_types_file (str): Absolute path of the domain types file
        """
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

    def read_in_known_PPIs(self):
        """
        Reads self.PPI_file and sets up protein pairs in tuples to store them in the self.known_PPIs set
        """
        with open(self.PPI_file, 'r') as file:
            lines= [line.strip() for line in file.readlines()] # PRS saved as .tsv
        for line in lines:
            tab= line.split('\t')
            PPI_instance= sorted(list([tab[0], tab[1]])) # a sorted protein pair as list
            self.known_PPIs.add(tuple(PPI_instance)) # PPI pair saved as tuple
        print(f'{len(self.known_PPIs)} PPIs read in.')

    def read_in_domain_matches(self,domain_matches_json_file): # This one reads in all domain matches, useful for DDI predictor
        """
        Reads a json file that contains domain matches in human proteins and sets up DomainMatch instances that are subsequently stored in the domain_matches_dict attribute of the protein instance. Protein instances assessed through self.proteins_dict.

        Args:
            domain_matches_json_file (str): Absolute path of the domain matches json file
        """
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
