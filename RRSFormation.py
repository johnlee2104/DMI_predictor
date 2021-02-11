import DMIDB

class RRSFormation:
    def __init__(self, RRS_version, dump_file= None, write_file= None):
        self.RRS_instances= []
        self.RRS_version= RRS_version
        if dump_file is None:
            self.dump_file= str(RRS_version) + '.pickle'
        else:
            self.dump_file= dump_file
        if write_file is None:
            self.write_file= str(RRS_version) + '.tsv'
        else:
            self.write_file= write_file

    def dump_RRSFormation_instance(self):
        from pickle import dump
        with open(self.dump_file, 'wb') as pickle_file: 
            dump(self, pickle_file)

    def write_out_RRS_instances(self): # take RRS_version as argument to construct the file name of the output file
        file= open(self.write_file, 'w')
        file.write('\t'.join((' ', 'Elm', 'interactorElm', 'ElmMatch', 'IUPredLong', 'IUPredShort', 'Anchor',
        'interactorDomain', 'Domain_ID1', 'DomainMatch1', 'Evalue1', 'Domain_ID2', 'DomainMatch2', 'Evalue2')))
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
        file.close()
