import os,sys
import subprocess
import math
import environment_variables

import pdb_extraction as pdb
import prepare_dock as d

import adf_engine as adf 
import gaussian_engine as g
import orca_engine as orca

from parser import Parser


             #         r_OA, e_OA, r_SA, e_SA, r_HD, e_HD, r_NA, e_NA,  r_N,  e_N, r_M_HD, e_M_HD
standard_set = {'V' : [ 2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0],
                'CR': [ 2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0],
                'FE': [ 2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0],
                'CO': [ 2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0],
                'NI': [ 2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0],
                'CU': [ 2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0],
                'Y' : [ 2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0],
                'MO': [ 2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0],
                'RU': [ 1.171,  5.261,  1.114,  15.632,  1.116,  0.572,  2.900,  20.035,  1.733,  5.488,  1.975,  3.459],
                'RH': [ 2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0],
                'PD': [ 2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0],
                'RE': [ 2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0],
                'OS': [ 2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0],
                'IR': [ 2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0],
                'PT': [ 2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0],
                'AU': [ 2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0,  2.0,  5.0]
            }

def docking(input_file, par=None, test_GA=False):
    
    if test_GA == False:
        par = Parser(input_file)
    else:
        pass

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

    if os.path.exists(par.name_ligand+'_c.xyz') == False:
        subprocess.call([os.environ['OBABEL']+f' -ixyz {input_dir}/'+par.xyz_file+' -oxyz '+par.name_ligand+'_c.xyz --canonical > '+par.name_ligand+'_c.xyz'],shell=True,  stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    if os.path.exists(par.name_ligand+'.mol2') == False:
        subprocess.call([os.environ['OBABEL']+' -ixyz '+par.name_ligand+'_c.xyz -omol2 '+par.name_ligand+'.mol2  > '+par.name_ligand+'.mol2'],shell=True,  stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    if os.path.exists(par.name_ligand+'.sdf') == False:
        subprocess.call([os.environ['OBABEL']+' -ixyz '+par.name_ligand+'_c.xyz -osdf '+par.name_ligand+'.sdf  > '+par.name_ligand+'.sdf'],shell=True,  stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    ###### Create pdb files ###### 
    os.system(f'cp {input_dir}/'+par.pdb_file+' .')

    if os.path.exists('clean_'+par.name_protein+'.pdb') == False:
        os.system(f"cp  {input_dir}/"+par.name_protein+".pdb .")
        pdb.protonate_pdb(par.pdb_file, par.pH)
        pdb.clean_protein_pdb(par.name_protein, par.pdb_file, par.clean_pdb)

    ###### Quantum Calculations ######
    os.chdir(f'{output_dir}')

    if os.path.isdir('QM') == False:
        os.mkdir('QM')
        os.chdir('QM')
    else:
        os.chdir('QM')

    if par.engine.lower() == 'adf':
        qm_dir, energy = adf.adf_engine(f'{output_dir}/file_prep/'+par.name_ligand+'_c.xyz', par, output_dir)

    if par.engine.lower() == 'gaussian':
        qm_dir, energy = g.gaussian_engine(f'{output_dir}/file_prep/'+par.name_ligand+'_c.xyz', par, output_dir)

    if par.engine.lower() == 'orca':
        qm_dir, energy = orca.orca_engine(f'{output_dir}/file_prep/'+par.name_ligand+'_c.xyz', par, output_dir)

    ##### AutoDock #####
    os.chdir(f'{output_dir}')

    if os.path.isdir('docking') == False:
        os.mkdir('docking')
        os.chdir('docking')
    else:
        os.chdir('docking')

    subprocess.call([f'cp {input_dir}/'+par.parameter_file+' .'], shell=True)

    subprocess.call([f'cp {output_dir}/file_prep/clean_'+par.name_protein+'.pdb .'], shell=True)

    subprocess.call([f'cp {output_dir}/file_prep/'+par.name_ligand+'.mol2 .'], shell=True)
    subprocess.call([f'cp {qm_dir}/CM5_charges .'], shell=True)

    subprocess.call([f'cp {output_dir}/file_prep/'+par.name_ligand+'_c.xyz .'], shell=True)

    if par.rmsd == True:
        if os.path.isfile('ref.xyz') == False:
            subprocess.call([os.environ['OBABEL']+" -ixyz "+par.name_ligand+"_c.xyz -oxyz ref.xyz -d > ref.xyz"], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    if par.dock_x and par.dock_y and par.dock_z != None:
        dock = d.users_coordinates(par.dock_x, par.dock_y, par.dock_z)
    else:
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

    if par.standard == True:
        parameter_set = standard_set.get(par.metal_symbol.upper())
        d.create_ligand_pdbqt_file(par.name_ligand)
        d.prepare_receptor(par.name_protein)
        d.docking_func(parameter_set, par.parameter_file, par.metal_symbol, par.name_ligand, par.name_protein, dock, box_size, par.num_poses, par.dock_algorithm, par.random_pos, par.ga_dock, par.sa_dock, energy)
    else:
        d.create_ligand_pdbqt_file(par.name_ligand)
        d.prepare_receptor(par.name_protein)
        d.docking_func(par.parameter_set, par.parameter_file, par.metal_symbol, par.name_ligand, par.name_protein, dock, box_size, par.num_poses, par.dock_algorithm, par.random_pos, par.ga_dock, par.sa_dock, energy)

    if par.rmsd == True:
        rmsd_list = []
        i = 1
        while os.path.exists(par.name_ligand+"_%i.pdbqt" % i):
            subprocess.call([os.environ['OBABEL']+" -ipdbqt "+par.name_ligand+"_{}.pdbqt".format(i)+" -oxyz "+par.name_ligand+"_{}.xyz".format(i)+" -d > "+par.name_ligand+"_{}.xyz".format(i)], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            rmsd_non_rotate = float(subprocess.getoutput([os.environ['PYTHON_3']+' '+os.environ['ROOT_DIR']+'/metal_dock/calculate_rmsd.py ref.xyz '+par.name_ligand+'_{}.xyz'.format(i)+' -nh --reorder --rotation none --translation none']))
            rmsd = rmsd_non_rotate

            rmsd_list.append("RMSD for Conformation %i = %.4f"% (i, rmsd))
            i += 1

        for i in range(0,len(rmsd_list)):
            print(rmsd_list[i])


    ##### results #####
    os.chdir(f'{output_dir}')

    if os.path.isdir('results') == False:
        os.mkdir('results')
        os.chdir('results')
    else:
        os.chdir('results')

    subprocess.call([f'mv {output_dir}/docking/'+par.name_ligand+'_*.pdbqt .'], shell=True)
    subprocess.call([f'mv {output_dir}/docking/clean_'+par.pdb_file+' .'], shell=True)
 
    print("\nDOCKING SUCCESFULLY COMPLETED")
    print("THE PRINTED POSES AND PROTEIN CAN BE FOUND IN THE RESULTS DIRECTORY")
    return