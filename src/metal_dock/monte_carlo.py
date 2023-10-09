
import os, shutil, glob, math, subprocess
import random
import numpy as np
import scipy as sc
from distutils.dir_util import copy_tree

from itertools import repeat
from . import prepare_dock as d
from multiprocessing import Pool
from scipy.stats import randint

def convertible(v):
    try:
        int(v)
        return True
    except (TypeError, ValueError):
        return False

def random_sample_continuous():
    return sc.random.uniform(low=0, high=7)

def dock_pool(n_prot, par, parameter_set, input_dir, tmp_dir):
    copy_tree(f'{input_dir}/data_set/protein_{n_prot}/output/docking', f'{os.getcwd()}/protein_{n_prot}/docking')
    os.chdir(f'{tmp_dir}/protein_{n_prot}/docking')
    output_temp_dir=f'{tmp_dir}/protein_{n_prot}/docking'

    # Obtain metal complex and protein names
    for files in glob.glob("*_c.xyz"):
        file_list = files.split('_c.xyz')
        name_ligand = file_list[0]

    for files in glob.glob("clean_*.pdb"):
        file_list = files.split('clean_')
        file_list = file_list[1].split('.')
        name_protein = file_list[0]

    ##### AutoDock ##### 
    dock = d.get_coordinates(f'{output_temp_dir}/ref.xyz',par.metal_symbol)
    
    if par.parameter_file == 'metal_dock.dat':
        shutil.copyfile(os.environ['ROOT_DIR']+'/metal_dock.dat', os.getcwd()+f'/metal_dock.dat')
        d.docking_func(par, parameter_set, name_ligand, name_protein, dock, par.box_size, energy=None)

    ##### Fitness function ######
    output =[]
    output.append(f'---------------------------------\n')
    output.append(f'PROTEIN {n_prot}\n')
    output.append(f'---------------------------------\n')
    rmsd_list = []
    i = 1
    while os.path.exists(f'{output_temp_dir}/{name_ligand}_{i}.pdbqt'):
        d.delete_hydrogen(f'{output_temp_dir}/{name_ligand}_{i}.pdbqt')
        subprocess.call([os.environ['OBABEL']+f' -ipdbqt {name_ligand}_{i}.pdbqt -oxyz {name_ligand}_{i}.xyz -d > {name_ligand}_{i}.xyz'], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        rmsd_non_rotate = float(subprocess.getoutput([os.environ['PYTHON_3']+' '+os.environ['ROOT_DIR']+'/metal_dock/calculate_rmsd.py ref.xyz '+name_ligand+'_{}.xyz'.format(i)+' -nh --reorder --rotation none --translation none']))
        rmsd = rmsd_non_rotate
        rmsd_list.append(rmsd)
        output.append(f"RMSD for Conformation {i} = {rmsd:.4f}\n")
        i += 1          

    avg_output = np.mean(rmsd_list)
    print(''.join(output))
    return avg_output


def calculate_rmsd(par, parameter_set, input_dir, tmp_dir, dir_list):
    # global parameter_average
    # parameter_average = []

    with Pool(processes=len(dir_list)) as pool:
        rmsd_avg_list = pool.starmap(dock_pool, zip(dir_list, repeat(par), repeat(parameter_set), repeat(input_dir), repeat(tmp_dir)))

    rmsd_avg = np.mean(np.array((rmsd_avg_list)))
    return rmsd_avg

def optimize_MC(input_file):
    input_dir = os.getcwd()
    output_dir = f'{input_dir}/output'
    
    par = input_file
    
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

    ###### Calculate Box Size #######
    if par.box_size != 0 and par.scale_factor == 0:
        par.box_size = par.box_size * 2.66 # Convert Å to grid points
    if int(par.box_size) == par.box_size:
        par.box_size =  int(par.box_size)
    else:
        par.box_size = math.ceil(par.box_size)
        print('SPACING BETWEEN GRID POINTS IS STANDARD SET TO 0.375 Å')
        print('BOX SIZE MUST BE INTEGER GRID POINTS WHICH WAS NOT FOUND')
        print('BOX SIZE SIDE ROUNDED UP AND SET TO {:.3f} Å\n'.format(par.box_size / 2.66))
    if par.box_size == 0 and par.scale_factor != 0:
        par.box_size = d.box_size_func(par.name_ligand+'_c.xyz', par.metal_symbol, 0.375, par.scale_factor)
    if par.box_size != 0 and par.scale_factor != 0:
        print("CANNOT SELECT BOXSIZE AND SCALE FACTOR - SET ONE VALUE TO 0")
        sys.exit()
    if par.box_size == 0 and par.scale_factor == 0:
        print("CANNOT SELECT BOXSIZE AND SCALE FACTOR - SET ONE VALUE GREATER THAN 0")
        sys.exit()

    if os.path.exists('parameter_history'):
        os.remove('parameter_history')

    with open('parameter_history', 'a') as f:
        f.write(f"PARAMETERS          :        e_NA        e_OA        e_SA        e_HD |       RMSD \n")

    if os.path.isdir(f'{output_dir}/tmp'):
        shutil.rmtree(f'{output_dir}/tmp', ignore_errors=True)
        os.mkdir(f'{output_dir}/tmp')
        os.chdir(f'{output_dir}/tmp')
    else:
        os.mkdir(f'{output_dir}/tmp')
        os.chdir(f'{output_dir}/tmp')

    tmp_dir=os.getcwd()

    os.chdir(f'{tmp_dir}')

    new_parameter_set = [2, 2, 2, 2]
    best_parameter_set = [2, 2, 2, 2]
    best_rmsd = calculate_rmsd(par, best_parameter_set, input_dir, tmp_dir, dir_list)

    i=0
    step=1
    while i < par.mc_steps:
        for idx in range(len(best_parameter_set)):

            print(f'######################################################################################\n')
            print(f'################################ MONTE CARLO STEP {step} ##################################\n')
            print(f'######################################################################################\n')

            random_eps = random_sample_continuous()
            new_parameter_set[idx] = random_eps

            rmsd = calculate_rmsd(par, new_parameter_set, input_dir, tmp_dir, dir_list)

            print(f'RMSD: {rmsd}')
            print(f'PARAMETER SET: {new_parameter_set}\n')

            print(f'BEST RMSD: {best_rmsd}')
            print(f'BEST PARAMETER SET {best_parameter_set}\n')
    
            if rmsd < best_rmsd:
                best_parameter_set[idx] = new_parameter_set[idx]
                best_rmsd = rmsd
                print('PARAMETER SET ACCEPTED')
                print(f'NEW BEST PARAMETER SET {best_parameter_set}\n')
                with open(f'{output_dir}/parameter_history', 'a') as f:
                    f.write('STEP {:>6}         :  '.format(step)+'  '.join(format(best_parameter_set[x], ">10.5f") for x in range(0,len(best_parameter_set)))+' | {:>10.5f}  \n'.format(best_rmsd))
            else:
                diff_RMSD = best_rmsd - rmsd 
                acceptance_probability = np.exp(diff_RMSD)
                random_number = random.uniform(0,1)

                if random_number < acceptance_probability:
                    best_parameter_set[idx] = new_parameter_set[idx]
                    best_rmsd = rmsd
                    print('ACCEPTED WITH HIGHER RMSD')
                    print(f'NEW PARAMETER SET {best_parameter_set}\n')
                    with open(f'{output_dir}/parameter_history', 'a') as f:
                        f.write('STEP {:>6}         :  '.format(step)+'  '.join(format(best_parameter_set[x], ">10.5f") for x in range(0,len(best_parameter_set)))+' | {:>10.5f}  \n'.format(best_rmsd))

                else:
                    new_parameter_set[idx] = best_parameter_set[idx]
                    print('PARAMETER SET DENIED')
            step+=1

        i+=1

    with open(f'{output_dir}/parameter_history','r') as fin:
        lines = [line.strip().split() for line in fin]
        rmsd =  [line[-1] for line in lines[1:]]
        min_idx = np.argmin(rmsd)+1
        line = lines[min_idx]
    
    shutil.rmtree(f'{tmp_dir}') 
    os.chdir(f'{output_dir}')
    print(f'#################################################################################################################################\n')
    print(f'BEST RMSD: {line[-1]}')
    print(f'BEST PARAMETERS: e_NA {line[3]} kcal/mol; e_OA {line[4]} kcal/mol; e_SA {line[5]} kcal/mol; e_HD {line[6]} kcal/mol')
