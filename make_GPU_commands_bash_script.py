# This script takes a directory e.g. run17 and finds all the fasta files in the directory to generate bash scripts that are used to deploy GPU for AF2 prediction, The number of GPU available needs to be specified in order to generate the number of bash scripts that evenly distribute the workload to the GPU available.

import os
import sys

def make_sh_scripts(dir_path, GPU_num):
    GPU_num= int(GPU_num)
    header = "#! /bin/bash"

    run_id = os.path.split(dir_path)[-1]
    fasta_list = []
    for file in os.listdir(dir_path):
        if os.path.splitext(file)[1] == '.fasta':
            fasta_list.append(file)

    for i in range(GPU_num):
        with open(f'{dir_path}/{run_id}_GPU{i + 1}_commands.sh', 'w') as f:
            f.write(f'{header}\n')
            for ind, fasta in enumerate(fasta_list):
                if (ind + 1) % GPU_num == i:
                    f.write('\n')
                    fasta = fasta
                    command = f"""SINGULARITYENV_CUDA_VISIBLE_DEVICES={i+1} time singularity run --nv /fsimb/common/singularity_tools/alphafold/v2.1.1r0/alphafold_2.1.1.simg \\\n--data_dir=/media/storage/alphafold_v211_databases \\\n--fasta_paths=/fsimb/groups/imb-luckgr/projects/dmi_predictor/DMI_AF2_PRS/{run_id}/{fasta} \\\n--max_template_date=2020-05-14 \\\n--model_preset=multimer \\\n--is_prokaryote_list=false \\\n--db_preset=full_dbs \\\n--output_dir=/fsimb/groups/imb-luckgr/projects/dmi_predictor/DMI_AF2_PRS/{run_id} \\\n--uniref90_database_path=/media/storage/alphafold_v211_databases/uniref90/uniref90.fasta \\\n--mgnify_database_path=/media/storage/alphafold_v211_databases/mgnify/mgy_clusters_2018_12.fa \\\n--template_mmcif_dir=/media/storage/alphafold_v211_databases/pdb_mmcif/mmcif_files \\\n--obsolete_pdbs_path=/media/storage/alphafold_v211_databases/pdb_mmcif/obsolete.dat \\\n--uniclust30_database_path=/media/storage/alphafold_v211_databases/uniclust30/uniclust30_2018_08/uniclust30_2018_08 \\\n--pdb_seqres_database_path=/media/storage/alphafold_v211_databases/pdb_seqres/pdb_seqres.txt \\\n--uniprot_database_path=/media/storage/alphafold_v211_databases/uniprot/uniprot.fasta \\\n--bfd_database_path=/media/storage/alphafold_v211_databases/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"""
                    f.write(f'{command}\n')

def main(dir_path, GPU_num):
    make_sh_scripts(dir_path, GPU_num)

if __name__ == '__main__':
    dir_path = sys.argv[1]
    GPU_num = sys.argv[2]
    main(dir_path, GPU_num)

# python3 make_GPU_commands_bash_script.py ./run17 7