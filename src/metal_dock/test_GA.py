import os, glob, subprocess

from parser import Parser
from docking import docking


def convertible(v):
    try:
        int(v)
        return True
    except (TypeError, ValueError):
        return False

def test_GA(input_file):
    par = Parser(input_file)

    par.rmsd = True

    input_dir = os.getcwd()
    output_dir = f'{input_dir}/output'

    os.chdir(f'{input_dir}/data_set')

    # Make list of the protein numbers to iterate over
    dir_list = os.listdir(os.getcwd())
    dir_list = [str(i).replace('protein_','') for i in dir_list]
    dir_list = [int(i) for i in dir_list if convertible(i)]
    dir_list = sorted(dir_list)

    os.chdir(f'{input_dir}')

    ###### Generate Output Dir #######
    if os.path.isdir('output') == False:
        os.mkdir('output')
        os.chdir('output')
    else:
        os.chdir('output')

    for n_prot in dir_list:
        subprocess.call(['cp','-r', f'{input_dir}/data_set/protein_{n_prot}', os.getcwd()], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        os.chdir(f'{output_dir}/protein_{n_prot}')

        # Obtain ligand and protein names
        for files in glob.glob("*.xyz"):
            par.xyz_file = files

            file_list = files.split('.xyz')
            par.name_ligand = file_list[0]


        for files in glob.glob("*.pdb"):
            par.pdb_file = files

            file_list = files.split('.pdb')
            par.name_protein = file_list[0]
        

        docking(input_file, par, test_GA=True)

        os.chdir(f'{output_dir}')