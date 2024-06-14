import os,sys, shutil
import subprocess
import math
import numpy as np

from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.io import read, write
from xtb.ase.calculator import XTB

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
    if os.path.isdir(output_dir) == False:
        os.mkdir(output_dir)
        os.chdir(output_dir)
    else:
        os.chdir(output_dir)

    if os.path.isdir(f'{output_dir}/file_prep') == False:
        os.mkdir(f'{output_dir}/file_prep')
        os.chdir(f'{output_dir}/file_prep')
    else:
        os.chdir(f'{output_dir}/file_prep')

    if os.path.exists(f'{par.name_ligand}_c.xyz') == False:
        xyz_file = os.path.join(input_dir, par.xyz_file)
        # center the molecule around the metal_symbol --> messes up the reference xyz file for rmsd 
        # xyz_file = d.center_molecule(input_dir, xyz_file, par.metal_symbol)
        subprocess.call([os.environ['OBABEL']+f' -ixyz {xyz_file} -oxyz {par.name_ligand}_c.xyz --canonical > {par.name_ligand}_c.xyz'],shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        subprocess.call([os.environ['OBABEL']+f' -ixyz {par.name_ligand}_c.xyz -opdb {par.name_ligand}_c.pdb > {par.name_ligand}_c.pdb'],shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    if os.path.exists(f'{par.name_ligand}.mol2') == False:
        subprocess.call([os.environ['OBABEL']+f' -ixyz {par.name_ligand}_c.xyz -omol2 {par.name_ligand}.mol2  > {par.name_ligand}.mol2'],shell=True,  stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    if os.path.exists(f'{par.name_ligand}.sdf') == False:
        subprocess.call([os.environ['OBABEL']+f' -ixyz {par.name_ligand}_c.xyz -osdf {par.name_ligand}.sdf  > {par.name_ligand}.sdf'],shell=True,  stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    ###### Create pdb files ###### 
    if os.path.exists(f'clean_{par.name_protein}.pdb') == False:
        input_pdb = os.path.join(input_dir, par.pdb_file)
        output_pdb = os.path.join(f'{output_dir}/file_prep', f'{par.name_protein}.pdb')
        
        shutil.copyfile(input_pdb, output_pdb)
        pdb.protonate_pdb(par.pdb_file, par.pH, par.clean_pdb)
        pdb.clean_protein_pdb(par.name_protein, par.clean_pdb)

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

    # if par.engine.lower() == 'adf':
    #     mol_in = os.path.join(output_dir, 'QM', 'geom_opt', 'plams_workdir', 'output.mol')
    # else:
    # mol_in = os.path.join(output_dir, 'QM', 'geom_opt', 'output.mol')
    # mol_out = os.path.join(os.getcwd(), f'{par.name_ligand}.mol')
    # shutil.copyfile(mol_in, mol_out)

    # f'{output_dir}/QM/geom_opt/output.mol'
    if par.geom_opt == True:
        subprocess.call([os.environ['OBABEL']+f' -imol {output_dir}/QM/geom_opt/output.mol -omol2 {par.name_ligand}.mol2  > {par.name_ligand}.mol2'],shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    else:
        subprocess.call([os.environ['OBABEL']+f' -imol {output_dir}/QM/single_point/output.mol -omol2 {par.name_ligand}.mol2  > {par.name_ligand}.mol2'],shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

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

    ##### results #####
    os.chdir(f'{output_dir}')

    if os.path.isdir('results') == False:
        os.mkdir('results')
        os.chdir('results')
    else:
        os.chdir('results')


    print('#==============================================================================#')
    print("ADDING AND OPTIMIZING HYDROGEN ATOMS TO THE METAL COMPLEX POSES")

    i = 1
    while os.path.exists(os.path.join(output_dir,'docking',f'{par.name_ligand}_{i}.pdbqt')):
        if os.path.isdir(f'{output_dir}/results/pose_{i}') == False:
            os.mkdir(f'{output_dir}/results/pose_{i}')
            os.chdir(f'{output_dir}/results/pose_{i}')
        else:
            os.chdir(f'{output_dir}/results/pose_{i}')

        pdqt_in = os.path.join(output_dir,'docking',f'{par.name_ligand}_{i}.pdbqt')
        pdqt_out = os.path.join(os.getcwd(), f'{par.name_ligand}_{i}.pdbqt')
        shutil.copyfile(pdqt_in, pdqt_out)

        # remove dummy atom if present 
        d.delete_dummy_atom(pdqt_out)

        subprocess.call([os.environ['OBABEL']+f" -ipdbqt {par.name_ligand}_{i}.pdbqt -oxyz {par.name_ligand}_{i}.xyz > {par.name_ligand}_{i}.xyz"], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # add hydrogens to the xyz file 
        atom_constraints = d.add_non_polar_hydrogens(f'{output_dir}/file_prep/{par.name_ligand}_c.xyz',
                                                        f'{output_dir}/file_prep/{par.name_ligand}_c.pdb',
                                                        f'{par.name_ligand}_{i}.xyz', 
                                                        f'{par.name_ligand}_{i}_H.xyz')
        
        print('\n#-------------------------------------------------#')
        print(f'Optimizing Hydrogen Atoms with GFN2-xTB of Pose {i}')
        atoms = read(f'{par.name_ligand}_{i}_H.xyz', format='xyz')
        atoms.calc = XTB(method='GFN2-xTB')
        atoms.set_constraint(FixAtoms(indices=atom_constraints))
        opt = QuasiNewton(atoms, trajectory=f'{par.name_ligand}_{i}_H.traj')
        opt.run(fmax=0.05)

        write(f'{par.name_ligand}_{i}_H.xyz', atoms, format='xyz')
        d.write_pose_to_pdb(f'{par.name_ligand}_{i}_H.xyz', f'{par.name_ligand}_{i}_H.pdb')
        # remove the xyz file and pdbqt file 
        os.remove(f'{par.name_ligand}_{i}.xyz')
        os.remove(f'{par.name_ligand}_{i}.pdbqt')
        i += 1

    # copy clean protein file to results directory
    os.chdir(f'{output_dir}/results')
    clean_in = os.path.join(output_dir,'docking',f'clean_{par.name_protein}.pdb')
    clean_out = os.path.join(os.getcwd(), f'clean_{par.name_protein}.pdb')
    shutil.copyfile(clean_in, clean_out)

    # copy the dlg file to results directory
    dlg_in = os.path.join(output_dir,'docking',f'{par.name_ligand}_clean_{par.name_protein}.dlg')
    dlg_out = os.path.join(os.getcwd(), f'docking_results.dlg')
    shutil.copyfile(dlg_in, dlg_out)

    binding_energy, ligand_efficiency = d.extract_dlg(dlg_out, par)
    print('\n#==============================================================================#')
    print('DOCKING RESULTS:')
    print('Ligand Efficiency = (binding energy) / (number of heavy atoms in metal complex)')
    print('Interacting Residues = residues within 4 Angstrom of the metal complex\n')

    for i in range(par.num_poses):
        print(f"Pose {i+1}:")
        print("-------------")
        pose_residues = d.extract_interacting_residues(f'{output_dir}/results/pose_{i+1}/{par.name_ligand}_{i+1}_H.xyz', f'clean_{par.name_protein}.pdb')
        print(f'Binding Energy: {binding_energy[i][1]:7.4f} kcal/mol')
        print(f'Ligand Efficiency: {ligand_efficiency[i]:7.4f} kcal/mol')
        print(f'Interacting Residues:')
        for residue in pose_residues:
            print(f'Residue: {residue[0]}, ID: {residue[1]:>3}')
        print('\n')

    if par.rmsd == True:
        print('\n#==============================================================================#')
        print("CALCULATING RMSD VALUES FOR EACH POSE WITH RESPECT TO THE STARTING XYZ FILE\n")

        rmsd_list = []
        avg_list = []
        rmsd_func = os.path.join(os.environ['ROOT_DIR'], 'metal_dock','calculate_rmsd.py')

        for pose in range(par.num_poses):
            os.chdir(f'{output_dir}/results/pose_{pose+1}')
            rmsd_non_rotate = float(subprocess.getoutput([os.environ['PYTHON_3']+f' {rmsd_func} {output_dir}/file_prep/{par.name_ligand}_c.xyz {output_dir}/results/pose_{pose+1}/{par.name_ligand}_{pose+1}_H.xyz -nh --reorder --rotation none --translation none']))
            rmsd = rmsd_non_rotate

            avg_list.append(rmsd)

            rmsd_list.append(f"RMSD for Conformation {pose+1:>3} = {rmsd:>8.4f}")

        for i in range(0,len(rmsd_list)):
            print(rmsd_list[i])

        avg_output = np.mean(avg_list)
        print(f'Average RMSD              = {avg_output:8.4f}')
        stdv_rmsd = np.std(avg_list)
        print(f'Standard Deviation RMSD   = {stdv_rmsd:8.4f}')
        var_rmsd = np.var(avg_list)
        print(f"Variance RMSD             = {var_rmsd:8.4f}\n")

    print('\n#==============================================================================#')
    print("METALDOCK SUCCESFULLY COMPLETED")
    print("THE PRINTED POSES AND PROTEIN CAN BE FOUND IN THE RESULTS DIRECTORY")
    print("EACH PDB FILE IN THE RESULTS/POSE_X DIRECTORY CAN BE VISUALIZED WITH E.G. PYMOL")


    # ##### results #####
    # os.chdir(f'{output_dir}')

    # if os.path.isdir('results') == False:
    #     os.mkdir('results')
    #     os.chdir('results')
    # else:
    #     os.chdir('results')

    # i = 1
    # while os.path.exists(os.path.join(output_dir,'docking',f'{par.name_ligand}_{i}.pdbqt')):
    #     pdqt_in = os.path.join(output_dir,'docking',f'{par.name_ligand}_{i}.pdbqt')
    #     pdqt_out = os.path.join(os.getcwd(), f'{par.name_ligand}_{i}.pdbqt')
    #     shutil.copyfile(pdqt_in, pdqt_out)
    #     if par.rmsd == True:
    #         xyz_in = os.path.join(output_dir,'docking',f'{par.name_ligand}_{i}.xyz')
    #         xyz_out = os.path.join(os.getcwd(), f'{par.name_ligand}_{i}.xyz')
    #         shutil.copyfile(xyz_in, xyz_out)
    #     i += 1
    
    # clean_in = os.path.join(output_dir,'docking',f'clean_{par.name_protein}.pdb')
    # clean_out = os.path.join(os.getcwd(), f'clean_{par.name_protein}.pdb')
    # shutil.copyfile(clean_in, clean_out)
