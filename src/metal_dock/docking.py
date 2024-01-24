import os,sys, shutil
import subprocess
import math
import numpy as np

from . import environment_variables

from . import pdb_extraction as pdb
from . import prepare_dock as d

from . import adf_engine as adf 
from . import gaussian_engine as g
from . import orca_engine as orca

def docking(input_file, par=None):

    par = input_file

    input_dir = os.getcwd()
    output_dir = input_dir+'/output'

    ###### Generate Output Dir #######
    if os.path.isdir('output') == False:
        os.mkdir('output')
        os.chdir('output')
    else:
        os.chdir('output')

    if os.path.isdir('file_prep') == False:
        os.mkdir('file_prep')
        os.chdir('file_prep')
    else:
        os.chdir('file_prep')

    if os.path.exists(f'{par.name_ligand}_c.xyz') == False:
        xyz_file = os.path.join(input_dir, par.xyz_file)
        subprocess.call([os.environ['OBABEL']+f' -ixyz {xyz_file} -oxyz {par.name_ligand}_c.xyz --canonical > {par.name_ligand}_c.xyz'],shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    if os.path.exists(f'{par.name_ligand}.mol2') == False:
        subprocess.call([os.environ['OBABEL']+f' -ixyz {par.name_ligand}_c.xyz -omol2 {par.name_ligand}.mol2  > {par.name_ligand}.mol2'],shell=True,  stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    if os.path.exists(f'{par.name_ligand}.sdf') == False:
        subprocess.call([os.environ['OBABEL']+f' -ixyz {par.name_ligand}_c.xyz -osdf {par.name_ligand}.sdf  > {par.name_ligand}.sdf'],shell=True,  stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    ###### Create pdb files ###### 

    if os.path.exists(f'clean_{par.name_protein}.pdb') == False:
        input_pdb = os.path.join(input_dir, par.pdb_file)
        output_pdb = os.path.join(os.getcwd(), par.pdb_file)
        shutil.copyfile(input_pdb, output_pdb)
        pdb.protonate_pdb(par.pdb_file, par.pH, par.clean_pdb)
        pdb.clean_protein_pdb(par.name_protein, par.pdb_file, par.clean_pdb)

    ###### Quantum Calculations ######
    os.chdir(output_dir)

    if os.path.isdir('QM') == False:
        os.mkdir('QM')
        os.chdir('QM')
    else:
        os.chdir('QM')

    xyz_file = os.path.join(output_dir,'file_prep', f'{par.name_ligand}_c.xyz')

    if par.engine.lower() == 'adf':
        qm_dir, energy = adf.adf_engine(xyz_file, par, output_dir)

    if par.engine.lower() == 'gaussian':
        qm_dir, energy = g.gaussian_engine(xyz_file, par, output_dir)

    if par.engine.lower() == 'orca':
        qm_dir, energy = orca.orca_engine(xyz_file, par, output_dir)

    ##### AutoDock #####
    os.chdir(output_dir)

    if os.path.isdir('docking') == False:
        os.mkdir('docking')
        os.chdir('docking')
    else:
        os.chdir('docking')
    
    if par.parameter_file == 'metal_dock.dat':
        in_file = os.path.join(os.environ['ROOT_DIR'], 'metal_dock.dat')
        out_file = os.path.join(os.getcwd(), 'metal_dock.dat')
        shutil.copyfile(in_file, out_file)
    else:
        in_file = os.path.join(input_dir, par.parameter_file)
        out_file = os.path.join(os.getcwd(), par.parameter_file)
        shutil.copyfile(in_file, out_file)

    clean_pdb_in = os.path.join(output_dir, 'file_prep', f'clean_{par.name_protein}.pdb')
    clean_pdb_out = os.path.join(os.getcwd(), f'clean_{par.name_protein}.pdb')
    shutil.copyfile(clean_pdb_in, clean_pdb_out)

    mol2_in = os.path.join(output_dir, 'file_prep', f'{par.name_ligand}.mol2')
    mol2_out = os.path.join(os.getcwd(), f'{par.name_ligand}.mol2')
    shutil.copyfile(mol2_in, mol2_out)

    cm5_in = os.path.join(qm_dir, 'CM5_charges')
    cm5_out = os.path.join(os.getcwd(), 'CM5_charges')
    shutil.copyfile(cm5_in, cm5_out)

    c_xyz_in = os.path.join(output_dir, 'file_prep', f'{par.name_ligand}_c.xyz')
    c_xyz_out = os.path.join(os.getcwd(), f'{par.name_ligand}_c.xyz')
    shutil.copyfile(c_xyz_in, c_xyz_out)
    
    if par.rmsd == True:
        if os.path.isfile('ref.xyz') == False:
            subprocess.call([os.environ['OBABEL']+f" -ixyz {par.name_ligand}_c.xyz -oxyz ref.xyz -d > ref.xyz"], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


    if par.dock_x and par.dock_y and par.dock_z != None:
        dock = d.users_coordinates(par.dock_x, par.dock_y, par.dock_z)
    else:
        dock = d.get_coordinates(f'{par.name_ligand}_c.xyz', par.metal_symbol)

    if par.box_size != 0 and par.scale_factor == 0:
        npts = par.box_size * 2.66 # Convert Å to grid points
        if int(npts) == npts:
            box_size =  int(npts)
        else:
            box_size = math.ceil(npts)
            print('SPACING BETWEEN GRID POINTS IS STANDARD SET TO 0.375 Å')
            print('BOX SIZE MUST BE INTEGER GRID POINTS WHICH WAS NOT FOUND')
            print('BOX SIZE SIDE ROUNDED UP AND SET TO {:.3f} Å\n'.format(box_size / 2.66))

    if par.box_size == 0 and par.scale_factor != 0:
        box_size = d.box_size_func(f'{par.name_ligand}_c.xyz', par.metal_symbol, 0.375, par.scale_factor)

    if par.box_size != 0 and par.scale_factor != 0:
        print("CANNOT SELECT BOXSIZE AND SCALE FACTOR - SET ONE VALUE TO 0")
        sys.exit()

    if par.box_size == 0 and par.scale_factor == 0:
        print("CANNOT SELECT BOXSIZE AND SCALE FACTOR - SET ONE VALUE GREATER THAN 0")
        sys.exit()

    d.create_ligand_pdbqt_file(par, par.name_ligand)
    if os.path.isfile(f'clean_{par.name_protein}.pdbqt') == False:
        d.prepare_receptor(par.name_protein)
    d.docking_func(par, par.name_ligand, par.name_protein, dock, box_size, energy)

    if par.rmsd == True:
        rmsd_list = []
        avg_list = []
        i = 1
        while os.path.exists(os.path.join(output_dir,'docking',f'{par.name_ligand}_{i}.pdbqt')):
            d.delete_hydrogen(os.path.join(output_dir,'docking',f'{par.name_ligand}_{i}.pdbqt'))
            subprocess.call([os.environ['OBABEL']+f" -ipdbqt {par.name_ligand}_{i}.pdbqt -oxyz {par.name_ligand}_{i}.xyz -d > {par.name_ligand}_{i}.xyz"], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            rmsd_func = os.path.join(os.environ['ROOT_DIR'], 'metal_dock','calculate_rmsd.py')

            rmsd_non_rotate = float(subprocess.getoutput([os.environ['PYTHON_3']+f' {rmsd_func} ref.xyz {par.name_ligand}_{i}.xyz -nh --reorder --rotation none --translation none']))
            rmsd = rmsd_non_rotate

            avg_list.append(rmsd)

            rmsd_list.append("RMSD for Conformation %i = %.4f"% (i, rmsd))
            i += 1

        for i in range(0,len(rmsd_list)):
            print(rmsd_list[i])

        avg_output = np.mean(avg_list)
        print(f'Average RMSD: {avg_output:.4f}')
        stdv_rmsd = np.std(avg_list)
        print(f'Standard Deviation RMSD: {stdv_rmsd:.4}')
        var_rmsd = np.var(avg_list)
        print(f"Variance RMSD: {var_rmsd:.4}\n")


    ##### results #####
    os.chdir(f'{output_dir}')

    if os.path.isdir('results') == False:
        os.mkdir('results')
        os.chdir('results')
    else:
        os.chdir('results')

    i = 1
    while os.path.exists(os.path.join(output_dir,'docking',f'{par.name_ligand}_{i}.pdbqt')):
        pdqt_in = os.path.join(output_dir,'docking',f'{par.name_ligand}_{i}.pdbqt')
        pdqt_out = os.path.join(os.getcwd(), f'{par.name_ligand}_{i}.pdbqt')
        shutil.copyfile(pdqt_in, pdqt_out)
        if par.rmsd == True:
            xyz_in = os.path.join(output_dir,'docking',f'{par.name_ligand}_{i}.xyz')
            xyz_out = os.path.join(os.getcwd(), f'{par.name_ligand}_{i}.xyz')
            shutil.copyfile(xyz_in, xyz_out)
        i += 1
    
    clean_in = os.path.join(output_dir,'docking',f'clean_{par.pdb_file}')
    clean_out = os.path.join(os.getcwd(), f'clean_{par.pdb_file}')
    shutil.copyfile(clean_in, clean_out)
 
    print("\nDOCKING SUCCESFULLY COMPLETED")
    print("THE PRINTED POSES AND PROTEIN CAN BE FOUND IN THE RESULTS DIRECTORY")
    return
