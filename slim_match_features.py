# Keeps all the function that calculates the SLiM match features and output it as files.
import sys
sys.path.append('/Users/johnlee/Coding/Python/DMI/')
from iupred2a import iupred2a_lib

out_path= '/Users/johnlee/Coding/Python/DMI/PRS_RRS_union_sequences_features/'

def calculate_iupred_long(input_file):
    with open(input_file, 'r') as f:
        lines= [line.strip() for line in f.readlines()]
    protein_id= lines[0]
    seq= lines[1]
    result_iupredlong= iupred2a_lib.iupred(seq, mode= 'long')
    result_iupredshort= iupred2a_lib.iupred(seq, mode= 'short')
    result_anchor= iupred2a_lib.anchor2(seq)
    with open(out_path + 'IUPred_long/' + protein_id + '_iupredlong.txt', 'w') as f:
        f.write(protein_id + '\n')
        for pos, residue in enumerate(seq):
            f.write(f'{pos+1}\t{residue}\t{result_iupredlong[pos]:.4f}\n')
        f.close()
    with open(out_path + 'IUPred_short/' + protein_id + '_iupredshort.txt', 'w') as f:
        f.write(protein_id + '\n')
        for pos, residue in enumerate(seq):
            f.write(f'{pos+1}\t{residue}\t{result_iupredlong[pos]:.4f}\n')
        f.close()
    with open(out_path + 'Anchor/' + protein_id + '_anchor.txt', 'w') as f:
        f.write(protein_id + '\n')
        for pos, residue in enumerate(seq):
            f.write(f'{pos+1}\t{residue}\t{result_iupredlong[pos]:.4f}\n')
        f.close()
    print(f'{protein_id} calculated')


if __name__ == '__main__':
    input_file= sys.argv[1]
    calculate_iupred_long(input_file)
