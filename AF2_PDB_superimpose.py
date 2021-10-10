# script to superimpose AlphaFold2 predicted structures with PDB structures
# author: Katja Luck

from pymol import cmd
from pymol import stored

AF2_path = '/Users/luck/IMB/projects/DMI_predictor/AF2/'
PRS_path = AF2_path + 'DMI_structure_PRS/'

file_to_run = 'AF2_DMI_PRS_annotation.tsv'

file1 = open(PRS_path + file_to_run,'r')
entries = file1.readlines()
file1.close()

target = open(PRS_path + 'output_' + file_to_run,'w')
target.write('dmi_name\tRMSD_domain\tRMSD_peptide\tpdb_id\tchain_domain\tchain_motif\tdomain_start\t' + \
             'domain_end\tmotif_start\tmotif_end\tnum_align_atoms_domain\talign_score_domain\tnum_align_resi_domain\n')

for line in entries[1:]:
    tab_list = str.split(line[:-1],'\t')
    DMIname = tab_list[0]
    pdb_id = tab_list[1]
    chain_domain = tab_list[2]
    chain_motif = tab_list[3]
    domain_start = tab_list[4]
    domain_end = tab_list[5]
    motif_start = tab_list[6]
    motif_end = tab_list[7]

    # load the 'real' structure and the predicted structure into pymol
    cmd.load(PRS_path + DMIname + '/' + pdb_id + '.pdb','real')
    cmd.load(PRS_path + DMIname + '/AF2_' + pdb_id + '.pdb','AF2')

    # copy the relevant residues from the real and AF2 structure into new objects for further manipulation
    cmd.create('real_wk', f'(real and chain {chain_domain} and resi {domain_start}-{domain_end}) or ' + \
                            f'(real and chain {chain_motif} and resi {motif_start}-{motif_end})')
    cmd.create('AF2_wk', '(AF2)')

    # alter chain names of peptide and domain of predicted structure to fit the chain names of the real structure
    # this is necessary for the rms_cur step later
    cmd.alter('AF2_wk and chain A', f'chain="{chain_domain}"')
    cmd.sort()
    # calculate the position of the motif residues to assign them to the peptide chain later
    stored.c_alpha = []
    cmd.iterate_state(1,f'AF2_wk and chain {chain_domain} and n. CA','stored.c_alpha.append([x,y,z])')
    len_residues = len(stored.c_alpha)
    print(len_residues)
    len_motif = int(motif_end) - int(motif_start)
    print(len_motif)
    # assign the motif residues in the AF structure to the same chain ID as in the real structure
    cmd.alter(f'AF2_wk and chain {chain_domain} and resi {len_residues-len_motif}-{len_residues}', f'chain="{chain_motif}"')
    cmd.sort()

    # alter the residue numbers in the AF2 chains to fit the residue numbering in the real structure, this is
    # needed for the rms_cur command later
    cmd.alter(f'AF2_wk and chain {chain_domain}',f'resi=(str(int(resi)+{int(domain_start)-1}))')
    cmd.sort()
    cmd.alter(f'AF2_wk and chain {chain_motif}',f'resi=(str(int(resi)-{(len_residues-len_motif-int(motif_start))}))')
    cmd.sort()

    # superimpose the real and AF2 structure on the domains and retrieve the RMSD and other superimposition quality measures
    domain_super_out = cmd.align(f'AF2_wk and chain {chain_domain} and resi {domain_start}-{domain_end}', \
                                 f'real_wk and chain {chain_domain} and resi {domain_start}-{domain_end}', \
                                 cycles=0, object='super_domain', mobile_state=1, target_state=1)
    RMSD_domain = domain_super_out[0]
    num_align_atoms_domain = domain_super_out[1]
    align_score_domain = domain_super_out[5]
    num_align_resi_domain = domain_super_out[6]

    # calculate the RMSD on the peptides without making a fit first, but based on the fit from superimposing the domains
    RMSD_peptide = cmd.rms_cur(f'AF2_wk and chain {chain_motif}', f'real_wk and chain {chain_motif}',\
                                mobile_state=1,target_state=1,matchmaker=0,cycles=0,object='peptide_super')

    # write out the results
    target.write('\t'.join([DMIname,str(RMSD_domain),str(RMSD_peptide),pdb_id,chain_domain,chain_motif,domain_start,domain_end,\
                motif_start,motif_end,str(num_align_atoms_domain),str(align_score_domain),str(num_align_resi_domain)]) + '\n')

    # visualize the superimposed structures in a nice way and save an image of it
    cmd.hide('all')
    cmd.show('cartoon','AF2_wk or real_wk')
    cmd.show('sticks',f'(AF2_wk and chain {chain_motif}) or (real_wk and chain {chain_motif})')
    cmd.hide('(hydro)')
    cmd.color('skyblue','real_wk')
    cmd.color('tv_orange','AF2_wk')
    cmd.color('atomic',f'AF2_wk and chain {chain_motif} and not elem C')
    cmd.color('atomic',f'real_wk and chain {chain_motif} and not elem C')
    cmd.png(PRS_path + DMIname + '/' + DMIname + '_AF2_' + pdb_id + '_superimpose.png',ray=1)

    # save this pymol session
    cmd.save(PRS_path + DMIname + '/' + DMIname + '_AF2_' + pdb_id + '_superimpose.pse')

target.close()
