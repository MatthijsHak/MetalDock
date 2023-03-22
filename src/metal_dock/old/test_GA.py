import os, glob,shutil, subprocess
import math

import numpy as np 
from distutils.dir_util import copy_tree
from parser_metal_dock import Parser

import prepare_dock as d


def convertible(v):
    try:
        int(v)
        return True
    except (TypeError, ValueError):
        return False

def test_GA(input_file):
    if isinstance(input_file, Parser) == True:
        par = input_file
        print('\nSTARTING TEST PROCEDURE WITH OBTAINED PARAMETERS')
    else:
        par = Parser(input_file)

    par.rmsd = True

    input_dir = os.getcwd()

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

    output_dir = f'{input_dir}/output'

    avg_all_prot_list = []

    for n_prot in dir_list:
        copy_tree(f'{input_dir}/data_set/protein_{n_prot}/output/docking', os.getcwd()+f'/protein_{n_prot}/docking')

        os.chdir(f'{output_dir}/protein_{n_prot}/docking')

        # Obtain ligand and protein names
        for files in glob.glob("*_c.xyz"):
            par.xyz_file = files

            file_list = files.split('_c.xyz')
            par.name_ligand = file_list[0]


        for files in glob.glob("clean_*.pdb"):
            par.pdb_file = files

            file_list = files.split('.pdb')
            par.name_protein = file_list[0]
            par.name_protein = par.name_protein[6:]
       
        if os.path.isfile('ref.xyz') == False:
            subprocess.call([os.environ['OBABEL']+" -ixyz "+par.name_ligand+"_c.xyz -oxyz ref.xyz -d > ref.xyz"], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        dock = d.get_coordinates(par.name_ligand+'_c.xyz', par.metal_symbol)

        if par.box_size != 0 and par.scale_factor == 0:
            npts = par.box_size * 2.66 # Convert Å to grid points
            if int(npts) == npts:
                box_size =  npts
            else:
                box_size = math.ceil(npts)
                print('SPACING BETWEEN GRID POINTS IS STANDARD SET TO 0.375 Å')
                print('BOX SIZE MUST BE INTEGER GRID POINTS WHICH WAS NOT FOUND')
                print('BOX SIZE SIDE ROUNDED UP AND SET TO {:.3f} Å\n'.format(box_size / 2.66))
        if par.box_size == 0 and par.scale_factor != 0:
            box_size = d.box_size_func(par.name_ligand+'_c.xyz', par.metal_symbol, 0.375, par.scale_factor)
        if par.box_size != 0 and par.scale_factor != 0:
            print("CANNOT SELECT BOXSIZE AND SCALE FACTOR - SET ONE VALUE TO 0")
            sys.exit()
        if par.box_size == 0 and par.scale_factor == 0:
            print("CANNOT SELECT BOXSIZE AND SCALE FACTOR - SET ONE VALUE GREATER THAN 0")
            sys.exit()

        if os.path.exists('clean_'+par.name_protein+'.pdb') == False:
            shutil.copyfile(f'{output_dir}/protein_{n_prot}/{par.name_protein}.pdb',os.getcwd()+f'/{par.name_protein}.pdb')
            pdb.protonate_pdb(par.pdb_file, par.pH)
            pdb.clean_protein_pdb(par.name_protein, par.pdb_file, par.clean_pdb)

        #d.create_ligand_pdbqt_file(par, par.name_ligand)
        if os.path.exists('clean_'+par.name_protein+'.pdbqt') == False:
            d.prepare_receptor(par.name_protein)
            
        d.docking_func(par, par.parameter_set, par.name_ligand, par.name_protein, dock, box_size)

        rmsd_list = []
        avg_list = []
        i = 1
        while os.path.exists(f'{output_dir}/protein_{n_prot}/docking/{par.name_ligand}_{i}.pdbqt'):
            d.delete_hydrogen(f'{output_dir}/protein_{n_prot}/docking/{par.name_ligand}_{i}.pdbqt')
            subprocess.call([os.environ['OBABEL']+f' -ipdbqt {par.name_ligand}_{i}.pdbqt -oxyz {par.name_ligand}_{i}.xyz -d > {par.name_ligand}_{i}.xyz'], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            rmsd_non_rotate = float(subprocess.getoutput([os.environ['PYTHON_3']+' '+os.environ['ROOT_DIR']+'/metal_dock/calculate_rmsd.py ref.xyz '+par.name_ligand+'_{}.xyz'.format(i)+' -nh --reorder --rotation none --translation none']))
            rmsd = rmsd_non_rotate

            avg_list.append(rmsd)

            rmsd_list.append("RMSD for Conformation %i = %.4f"% (i, rmsd))
            i += 1

        print(f'-------------------------------------------     PROTEIN {n_prot}      --------------------------------------------')
        for i in range(0,len(rmsd_list)):
            print(rmsd_list[i])       
            
        avg_list = np.array(avg_list)

        avg_output = np.mean(avg_list)
        print(f'Average RMSD protein_{n_prot}: {avg_output:.4f}')
        min_rmsd = np.min(avg_list)
        print(f"Lowest RMSD protein_{n_prot} {min_rmsd:.4f}")
        stdv_rmsd = np.std(avg_list)
        print(f'Standard Deviation RMSD: {stdv_rmsd:.4}')
        var_rmsd = np.var(avg_list)
        print(f"Variance RMSD: {var_rmsd:.4}\n")
        
        avg_all_prot_list.append(avg_output)

        os.chdir(f'{output_dir}')

        
    avg = np.mean(np.array((avg_all_prot_list)))

    print('\n')
    print('#' * 100)
    print("\nTEST GA SUCCESFULLY COMPLETED")
    print(f"TEST RMSD AVERAGE OF ALL PROTEINS IS: {avg}")
