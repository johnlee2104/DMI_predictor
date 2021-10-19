# script that receives an input file with a list of DMI types for
# which it should extract the domain and motif sequences from the pdb fasta
# file and create an input file for AF2 for each DMI type and PDB structure
# the files should be created in a run directory with the given ID
#

# input:
# run ID

import argparse
import os

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-r',help='Indicate the run number.')

    args = parser.parse_args()

    path = '/Volumes/imb-luckgr/projects/dmi_predictor/DMI_AF2_PRS/'

    file1 = open(path + 'run' + args.r + '/run' + args.r + '_DMI_type_list.txt','r')
    entries = file1.readlines()
    file1.close()

    linker = 'GGSGGSGGSGGSGGSGGSGGSGGSGGSGGS'

    for line in entries:
        DMI_type = line[:-1]
        fasta_path = path + 'DMI_types/' + DMI_type + '/'
        files = os.listdir(fasta_path)
        for fname in files:
            if fname[-6:] == '.fasta':
                pdb_code = fname[:-6]
                target = open(path + 'run' + args.r + '/run' + args.r + '_' + DMI_type + '_' + pdb_code + '.fasta','w')
                target.write('>run' + args.r + '_' + DMI_type + '_' + pdb_code + '\n')
                file2 = open(fasta_path + fname,'r')
                entries2 = file2.readlines()
                file2.close()
                seqs = []
                seq = ''
                for line2 in entries2[1:]:
                    if line2[0] == '>':
                        seqs.append(seq)
                        seq = ''
                    else:
                        seq = seq + line2[:-1]
                seqs.append(seq)
                print(DMI_type,pdb_code,seqs)
                AF2_seq = (seqs[0] + linker + seqs[1]) if len(seqs[0]) > len(seqs[1]) else (seqs[1] + linker + seqs[0])
                target.write(AF2_seq + '\n')
                target.close()
