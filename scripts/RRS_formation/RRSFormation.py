# This script contains the basal class for the formation of different versions of random reference set (RRS). The different versions of RRS are used to train random forest model and their performances are tested and compared to check potential biases in the model.
# Author: Chop Yan Lee

import DMIDB

class RRSFormation:
    """
    Represents a basal and generic class to form different versions of random reference set

    Attributes:
        RRS_instances (list of DMIMatch): Random DMI matches saved as DMIMatch instances in a list
        RRS_version (str): a version name given to RRS, e.g. RRSv1_1_20210427 as RRS version 1, triplicate 1 and date of generating the RRS
        dump_file (str): Absolute path that includes the name of a pickle file to be created, where the RRSFormation instance will be dumped into
        write_file (str): Name to be given to the file into which RRS instances are written out
    """

    def __init__(self, RRS_version, dump_file= None, write_file= None):
        """
        Instantiate RRSFormation

        Args:
            RRS_version (str): a version name given to RRS, e.g. RRSv1_1_20210427 as RRS version 1, triplicate 1 and date of generating the RRS
            dump_file (str): Absolute path that includes the name of a pickle file to be created, where the RRSFormation instance will be dumped into
            write_file (str): Name to be given to the file into which RRS instances are written out
        """
        self.RRS_instances= []
        self.RRS_version= RRS_version
        if dump_file == None:
            self.dump_file= str(RRS_version) + '.pickle'
        else:
            self.dump_file= dump_file
        if write_file == None:
            self.write_file= str(RRS_version) + '.tsv'
        else:
            self.write_file= write_file

    def dump_RRSFormation_instance(self):
        """
        Dump an instance of RRSFormation into a pickle file
        """
        from pickle import dump
        with open(self.dump_file, 'wb') as pickle_file:
            dump(self, pickle_file)
        print(f'{self.RRS_version} is pickled as {self.dump_file}.')

    def write_out_RRS_instances(self, InterfaceHandling): # take RRS_version as argument to construct the file name of the output file
        """
        Write out the DMIMatch instances stored in self.RRS_instances as random reference datapoints for model training and testing

        Args:
            InterfaceHandling (DMIDB.InterfaceHandling): An instance of the InterfaceHandling class
        """
        file= open(self.write_file, 'w')
        file.write('\t'.join((' ', 'Accession', 'Elm', 'Regex', 'Pattern', 'Probability', 'interactorElm', 'ElmMatch', 'interactorDomain', 'DomainID1', 'DomainMatch1', 'DomainMatchEvalue1', 'DomainFreqbyProtein1', 'DomainFreqinProteome1', 'DomainID2', 'DomainMatch2', 'DomainMatchEvalue2', 'DomainFreqbyProtein2', 'DomainFreqinProteome2')))
        file.write('\n')
        for i , inst in enumerate(self.RRS_instances):
            file.write('\t'.join((str(i+1), inst.slim_match.slim_type_inst.slim_id, inst.slim_match.slim_type_inst.name, inst.slim_match.slim_type_inst.regex, inst.slim_match.pattern, inst.slim_match.slim_type_inst.probability, inst.slim_protein, '-'.join([str(inst.slim_match.start) , str(inst.slim_match.end)]), inst.domain_protein)))
            for domain_match_list in inst.domain_interface_match.domain_matches:
                domain_id= domain_match_list[0].domain_id
                start= [domain_match.start for domain_match in domain_match_list]
                end= [domain_match.end for domain_match in domain_match_list]
                evalue= [domain_match.evalue for domain_match in domain_match_list]
                match= '|'.join([f"{i[0]}-{i[1]}" for i in zip(start, end)])
                evalues= '|'.join([str(ev) for ev in evalue])
                domainfreqbyprotein= InterfaceHandling.domain_types_dict[domain_id].DomainFreqbyProtein
                domainfreqinproteome= InterfaceHandling.domain_types_dict[domain_id].DomainFreqinProteome
                file.write('\t')
                file.write('\t'.join((domain_id, match, evalues, domainfreqbyprotein, domainfreqinproteome)))
            file.write('\n')
        file.close()
        print(f'{self.RRS_version} file saved as {self.write_file}.')
