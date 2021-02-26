# Keeps all the function that calculates the SLiM match features and output it as files.
import sys, json, os, glob
import numpy as np
sys.path.append('/Users/chopyanlee/Coding/Python/DMI/')
from iupred2a import iupred2a_lib

def calculate_iupred_scores(input_file, out_path):
    for file_name in glob.glob(in_path + '*'):
        with open(file_name, 'r') as f:
            lines= [line.strip() for line in f.readlines()]
        protein_id= lines[0][1:]
        seq= lines[1]
        result_iupredlong= iupred2a_lib.iupred(seq, mode= 'long')
        result_iupredshort= iupred2a_lib.iupred(seq, mode= 'short')
        result_anchor= iupred2a_lib.anchor2(seq)
        with open(out_path + '/IUPred_long/' + protein_id + '_iupredlong.txt', 'w') as f:
            f.write(protein_id + '\n')
            for pos, residue in enumerate(seq):
                f.write(f'{pos+1}\t{residue}\t{result_iupredlong[pos]}\n')
            f.close()
        with open(out_path + '/IUPred_short/' + protein_id + '_iupredshort.txt', 'w') as f:
            f.write(protein_id + '\n')
            for pos, residue in enumerate(seq):
                f.write(f'{pos+1}\t{residue}\t{result_iupredshort[pos]}\n')
            f.close()
        with open(out_path + '/Anchor/' + protein_id + '_anchor.txt', 'w') as f:
            f.write(protein_id + '\n')
            for pos, residue in enumerate(seq):
                f.write(f'{pos+1}\t{residue}\t{result_anchor[pos]}\n')
            f.close()
        print(f'IUPred scores and Anchor scores of {protein_id} calculated.')

def calculate_domain_overlap_score(in_path, out_path):
    for file_name in glob.glob(in_path + '*'):
        with open(file_name, 'r') as f:
            lines= [line.strip() for line in f.readlines()]
        protein_id= lines[0][1:]
        seq= lines[1]
        domain_overlap_scores= np.array([0] * len(seq))
        smart_pfam_domain_matches= [smart_domain_matches, pfam_domain_matches]
        for data in smart_pfam_domain_matches:
            for result in data['results']:
                if result['metadata']['accession']== protein_id:
                    for domain_match_id in result['entry_subset']:
                        for domain_match in domain_match_id['entry_protein_locations']:
                            for fragment in domain_match['fragments']:
                                domain_overlap_scores[int(fragment['start'])-1:int(fragment['end'])]= 1
        if protein_id in isoform_domain_matches:
            for domain_match in isoform_domain_matches[protein_id]:
                domain_overlap_scores[int(domain_match[1])-1:int(domain_match[2])]= 1
        with open(out_path + '/Domain_overlap/' + protein_id + '_domain_overlap.txt', 'w') as f:
            f.write(protein_id + '\n')
            for pos, residue in enumerate(seq):
                f.write(f'{pos+1}\t{residue}\t{domain_overlap_scores[pos]}\n')
            f.close()
        print(f'Domain overlap scores of {protein_id} calculated.')



if __name__ == '__main__':
    in_path= sys.argv[1]
    out_path= sys.argv[2]
    with open('/Users/chopyanlee/Coding/Python/DMI/domain_stuffs/interpro_9606_smart_matches_20210122.json', 'r') as f:
        smart_domain_matches= json.load(f)
    print('smart_json loaded')
    with open('/Users/chopyanlee/Coding/Python/DMI/domain_stuffs/interpro_9606_pfam_matches_20210122.json', 'r') as f:
        pfam_domain_matches= json.load(f)
    print('pfam_json loaded')
    with open('/Users/chopyanlee/Coding/Python/DMI/domain_stuffs/smart_domain_matches_for_alt_isoforms.txt', 'r') as f:
        lines= [line.strip() for line in f.readlines()]
    print('isoform_matches loaded')
    isoform_domain_matches= {}
    for line in lines:
        if line != '':
            if line[0] != '-':
                tab= [tab.strip() for tab in line.split('=')]
                if tab[0] == 'USER_PROTEIN_ID':
                    protein_id= tab[1]
                    if protein_id not in isoform_domain_matches:
                        isoform_domain_matches[protein_id]= []
                elif tab[0]== 'DOMAIN':
                    match= []
                    match.append(tab[1])
                elif tab[0]== 'START':
                    match.append(tab[1])
                elif tab[0]== 'END':
                    match.append(tab[1])
                    isoform_domain_matches[protein_id].append(match)
    for k, v in isoform_domain_matches.items():
        for ele in v:
            if 'coiled_coil_region' in ele:
                v.remove(ele)
    for line in lines:
        if line[:2] != '--':
            tab= line.split('=')
    calculate_iupred_scores(in_path, out_path)
    calculate_domain_overlap_score(in_path, out_path)
