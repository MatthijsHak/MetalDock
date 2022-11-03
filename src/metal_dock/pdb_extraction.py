import os,sys, shutil
import subprocess

def protonate_pdb(pdb_file, pH):
    subprocess.call(os.environ['PDB2PQR']+' --pdb-output pdb_prot.pdb --pH '+str(pH)+' --drop-water '+pdb_file+' pdb_prot.pdb', shell=True,  stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return

# Get clean.pdb
def clean_protein_pdb(name_protein, pdb_file, clean_pdb=True):
    if clean_pdb == True:
        # grep delete
        with open('pdb_prot.pdb', 'r') as fin:
            with open(f'clean_{pdb_file}', 'w') as fout:
                for line in fin:
                    if 'HETATM'  not in line:
                        fout.write(line)

    return