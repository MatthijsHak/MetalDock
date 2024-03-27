import os,sys, shutil
import subprocess

def protonate_pdb(pdb_file, pH, clean_pdb=True):
    if clean_pdb == True:
        subprocess.call(os.environ['PDB2PQR']+f' --noopt --pdb-output pdb_prot.pdb --with-ph {str(pH)} --drop-water {pdb_file} pdb_prot.pdb', shell=True,  stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    else:
        subprocess.call(os.environ['PDB2PQR']+f' --noopt --pdb-output pdb_prot.pdb --with-ph {str(pH)} {pdb_file} pdb_prot.pdb', shell=True,  stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return

# Get clean.pdb
def clean_protein_pdb(name_protein, clean_pdb=True):
    if clean_pdb == True:
        # grep delete
        with open('pdb_prot.pdb', 'r') as fin:
            with open(f'clean_{name_protein}', 'w') as fout:
                for line in fin:
                    if 'HETATM'  not in line:
                        fout.write(line)
    else:
        shutil.move('pdb_prot.pdb', f'clean_{name_protein}')

    return
