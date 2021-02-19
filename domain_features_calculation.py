import sys, json

""" This script takes a domain type files and returns a new file that contains annotation of the frequency of the
 domain in the human proteome."""

def calculate_domain_frequency(domain_type_file):
    with open(domain_type_file, 'r') as f:
        lines= [line.strip() for line in f.readlines()]
    total_human_protein= 20396
    if 'smart' in domain_type_file:
        data= smart_domain_matches
    elif 'pfam' in domain_type_file:
        data= pfam_domain_matches
    domain_match_in_protein= {} # No. of protein with at least 1 domain match: 2 repeat matches in 1 protein = 1
    total_domain_match_in_proteome= {} # Total number of domain match in the human proteome: 2 repeat matches in 1 protein = 2
    for line in lines[2:]:
        tab= line.split('\t')
        domain_id= tab[1]
        for result in data['results']:
            for domain_match_id in result['entry_subset']:
                if domain_id in domain_match_id.values():
                    domain_match_in_protein[domain_id]= domain_match_in_protein.get(domain_id, 0) + 1
                    total_domain_match_in_proteome[domain_id]= total_domain_match_in_proteome.get(domain_id, 0) + len(domain_match_id['entry_protein_locations'])
                    print(f'Frequency of {domain_id} calculated.')
    with open(domain_type_file[:-4] + '_with_frequency.txt', 'w') as f:
        f.write('\t'.join((lines[0], 'ProteinMatch', 'ProteomeMatch')))
        f.write('\n')
        f.write(lines[1])
        f.write('\n')
        for line in lines[2:]:
            tab= line.split('\t')
            domain_id= tab[1]
            if domain_id in domain_match_in_protein:
                f.write('\t'.join((line, str(domain_match_in_protein[domain_id]/total_human_protein),
                str(total_domain_match_in_proteome[domain_id]/total_human_protein))))
                f.write('\n')
            else:
                f.write('\t'.join((line, 'NA', 'NA'))) # Some domains are found only in certain taxon, and not human proteome
                f.write('\n')
        print(f'New file saved as {domain_type_file[:-4]}_with_frequency.txt')

if __name__ == '__main__':
    domain_type_file= sys.argv[1]
    with open('/Users/chopyanlee/Coding/Python/DMI/domain_stuffs/interpro_9606_smart_matches_20210122.json', 'r') as f:
        smart_domain_matches= json.load(f)
    if any(smart_domain_matches):
        print('smart_json loaded')
    with open('/Users/chopyanlee/Coding/Python/DMI/domain_stuffs/interpro_9606_pfam_matches_20210122.json', 'r') as f:
        pfam_domain_matches= json.load(f)
    if any(smart_domain_matches):
        print('smart_json loaded')
    calculate_domain_frequency(domain_type_file)
