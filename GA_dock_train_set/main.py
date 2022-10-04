#!/usr/bin/env python3 

import argparse
import glob,sys,os
import subprocess
import numpy as np
import pygad
import uuid
import shutil
import random
import statistics
import multiprocessing as mp 

from scipy.stats import rankdata

import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter

from multiprocessing import Pool
from distutils.dir_util import copy_tree


from openbabel import pybel as py
from openbabel import openbabel as ob

import docking as d
import input_variables as iv
import environment_variables

ob.obErrorLog.SetOutputLevel(0)

def convertible(v):
    try:
        int(v)
        return True
    except (TypeError, ValueError):
        return False

def is_float(string):
    try:
        float(string)
        return True
    except ValueError:
        return False

def flatten(l):
    return [item for sublist in l for item in sublist]

def mutation_func(offspring, ga_instance):
    '''Mutation function
    '''
    mutation_step = 100 
    mutation_positive = [x * 0.1 for x in range(1, mutation_step)]
    mutation_negative = [x * 0.001 for x in range(1, mutation_step)]

    #Adaptive Mutation
    mutation_probability_list = [ x  for x in np.linspace(0, 1, num=iv.var.sol_per_pop)]

    rank_list = rankdata(ga_instance.last_generation_fitness)

    offspring_length = iv.var.sol_per_pop - iv.var.keep_parents

    for i in range(0,offspring_length):
        mutation_probability = mutation_probability_list[int(rank_list[i])-1]

        for chromosome_idx in range(0,len(gene_space)):
            random_number_1 = random.uniform(0,1)
            random_number_2 = random.uniform(0,1)
            random_int = random.randrange(0, len(mutation_positive))

            if random_number_1 < mutation_probability:
                if random_number_2 <= 0.5:
                    mutated_gene = offspring[i][chromosome_idx] * mutation_positive[random_int]
                else:
                    mutated_gene = offspring[i][chromosome_idx] * mutation_negative[random_int]

                if gene_space[chromosome_idx]['low'] <= mutated_gene <= gene_space[chromosome_idx]['high']:
                    pass
                else:
                    mutated_gene = offspring[i][chromosome_idx]
            
            else:
                offspring[i] = offspring[i]

    return offspring

def box_size_func(sdf_file, spacing, scale_factor):
    # Open SDF file
    sdf = glob.glob(sdf_file)
    sdf = ''.join(str(x) for x in sdf)
    sdf = open(sdf, 'r')

    # Extract coordinates 
    lines = [line.split() for line in sdf]
    del lines[:4]

    for j in range(0,len(lines)):
        del lines[j][4:]
        string = lines[j][3]

        if is_float(string) == True:
            del lines[j:]
            break

    coordinates = []
    elements = []
    x_axis = []
    y_axis = []
    z_axis = []

    for k in range(0,len(lines)):
        if lines[k][3] == ''+iv.var.metal_symbol+'':
            metal = lines[k][:3]

        coordinates.append(lines[k][:3])
        elements.append(lines[k][3])

        x_axis.append(float(lines[k][0]))
        y_axis.append(float(lines[k][1]))
        z_axis.append(float(lines[k][2]))
        
    # Shift axis to centre at metal
    metal = [float(i) for i in metal]
    metal = np.array(metal)
    
    x_max = np.max(x_axis-metal[0])
    x_min = np.min(x_axis-metal[0])
    
    y_max = np.max(y_axis-metal[1])
    y_min = np.min(y_axis-metal[1])

    z_max = np.max(z_axis-metal[2])
    z_min = np.min(z_axis-metal[2])


    x_axis = np.abs(np.max(x_axis-metal[0]) - np.min(x_axis-metal[0]))*scale_factor
    y_axis = np.abs(np.max(y_axis-metal[1]) - np.min(y_axis-metal[1]))*scale_factor
    z_axis = np.abs(np.max(z_axis-metal[2]) - np.min(z_axis-metal[2]))*scale_factor

    if x_axis > 20:
        x_axis = 20
    if y_axis > 20:
        y_axis = 20
    if z_axis > 20:
        z_axis = 20

    x_npts = (round(x_axis / spacing)) & (-2)
    y_npts = (round(y_axis / spacing)) & (-2)
    z_npts = (round(z_axis / spacing)) & (-2)

    max_side = max([x_npts,y_npts,z_npts])
    print('Box size is {} {} {}'.format(x_npts*spacing,y_npts*spacing,z_npts*spacing))
    # Box dimensions in npts
    npts = [max_side, max_side, max_side]

    return npts

def docking_centre(coordinate_file):
    dock_site = open(coordinate_file,'r')
    coord = [line.split() for line in dock_site]

    if is_float(coord[0][0]) == True:
        dock = [coord[0][0], coord[0][1], coord[0][2]]
    else:
        dock = [coord[0][1], coord[0][2], coord[0][3]]

    return dock


def docking_func(solution, name_ligand, name_protein, energy, dock, npts, solution_idx):
    os.system('cp '+os.environ['WORKING_DIR']+'/'+iv.var.parameter_file+' .')
    
    # insert solution for R and epsilon for H-bond
    os.system(r'''awk '{ if ($2 == "'''+iv.var.metal_cap+'''" || $2 == "'''+iv.var.metal_symbol+'''") ($7 = '''+str(solution[10])+''') && ($8 = '''+str(solution[11])+'''); print $0}' '''+iv.var.parameter_file+''' > file_1''')
    os.system(r'''awk '{ if ($2 == "'''+iv.var.metal_cap+'''" || $2 == "'''+iv.var.metal_symbol+r'''") printf"%-8s %-3s %7s %8s %8s %9s %4s %4s %2s %3s %3s %2s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12; else print $0}' file_1 > '''+iv.var.parameter_file)
   
    #create_gpf():
    os.system(os.environ['PYTHON_2']+" "+os.environ['MGLTOOLS']+"/prepare_gpf4.py -l "+name_ligand+".pdbqt  -r clean_"+name_protein+".pdbqt -p parameter_file="+iv.var.parameter_file+" -p npts='{},{},{}'".format(npts[0],npts[1],npts[2])+" -p gridcenter='{:.4},{:.4},{:.4}' ".format(dock[0],dock[1],dock[2])+" >  /dev/null 2>&1")
    gpf = open('clean_'+name_protein+'.gpf', 'a')
    gpf.write('nbp_r_eps {:>.4f}  {:>.4f}'.format(solution[0],solution[1])+'    12 6 OA '+iv.var.metal_symbol+'\n')
    gpf.write('nbp_r_eps {:>.4f}  {:>.4f}'.format(solution[2],solution[3])+'    12 6 SA '+iv.var.metal_symbol+'\n')
    gpf.write('nbp_r_eps {:>.4f}  {:>.4f}'.format(solution[4],solution[5])+'    12 6 HD '+iv.var.metal_symbol+'\n')
    gpf.write('nbp_r_eps {:>.4f}  {:>.4f}'.format(solution[6],solution[7])+'    12 6 NA '+iv.var.metal_symbol+'\n')
    gpf.write('nbp_r_eps {:>.4f}  {:>.4f}'.format(solution[8],solution[9])+'    12 6  N '+iv.var.metal_symbol+'\n')
    gpf.close()

    #autogrid()
    os.system(os.environ['AUTODOCK']+'/autogrid4 -p clean_'+name_protein+'.gpf  >  /dev/null 2>&1')

    #create_dpf()
    d.write_dpf_file('clean_'+name_protein+'.gpf', name_ligand, 'clean_'+name_protein, iv.var.parameter_file, energy, random_pos=iv.var.random_position, SA=iv.var.docking_simulated_annealing, GA=iv.var.docking_genetic_algorithm)

    #autodock()
    os.system(os.environ['AUTODOCK']+'/autodock4 -p '+name_ligand+'_clean_'+name_protein+'.dpf  >  /dev/null 2>&1')

    #write_all_conformations()
    os.system(os.environ['PYTHON_2']+" "+os.environ['MGLTOOLS']+"/write_conformations_from_dlg.py -d "+name_ligand+"_clean_"+name_protein+".dlg")
    print('Parent {}: succesful!'.format(solution_idx))
    return

def rmsd_func(name_ligand, n_prot):
    rmsd_avg = []
    rmsd_list = []
    avg_list = []
    min_list = []
    rmsd_print_list = []
                
    output = [f"-------------------------------------------     PROTEIN {n_prot}      --------------------------------------------\n"]
              
    i = 1
    while os.path.exists(name_ligand+"_%i.pdbqt" % i):
        subprocess.call(os.environ['OBABEL']+" -ipdbqt "+name_ligand+"_{}.pdbqt".format(i)+" -oxyz "+name_ligand+"_{}.xyz".format(i)+" -d > "+name_ligand+"_{}.xyz".format(i), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) 

        rmsd_non_rotate = float(subprocess.getoutput([os.environ['PYTHON_3']+' '+os.environ['LIB_DIR']+'/calculate_rmsd.py ref.xyz '+name_ligand+'_{}.xyz'.format(i)+' --reorder --rotation none --translation none']))
        rmsd = rmsd_non_rotate

        rmsd_print_list.append("RMSD for Conformation %i = %.4f\n"% (i, rmsd))
        rmsd_list.append(rmsd)
        i += 1

    for j in range(0,len(rmsd_print_list)):
        output.append(rmsd_print_list[j])

        if generation == 0:
            first_gen = open(os.environ['WORKING_DIR']+'/first_gen/all_conf_first_gen','a')
            first_gen.write(rmsd_print_list[j])

            protein = open(os.environ['WORKING_DIR']+f'/first_gen/protein_{n_prot}','a')
            protein.write(rmsd_print_list[j])

        if generation == iv.var.num_generations+1:   
            last_gen = open(os.environ['WORKING_DIR']+'/last_gen/all_conf_last_gen','a')
            last_gen.write(rmsd_print_list[j])
            
            protein = open(os.environ['WORKING_DIR']+f'/last_gen/protein_{n_prot}','a')
            protein.write(rmsd_print_list[j])

    avg_output = np.mean(rmsd_list)
    avg_list.append(avg_output)
    output.append(f"Average RMSD: {avg_output:.4}\n")

    minimum_rmsd = min(rmsd_list)
    min_list.append(minimum_rmsd)
    output.append(f"Lowest RMSD: {minimum_rmsd:.4}\n")

    stdv_rmsd = np.std(rmsd_list)
    output.append(f"Standard Deviation RMSD: {stdv_rmsd:.4}\n")

    var_rmsd = np.var(rmsd_list)
    output.append(f"Variance RMSD: {var_rmsd:.4}\n")
    output.append("-----------------------------------------------------------------------------------------------------------\n")

    return avg_list, min_list, print(''.join(output))
    


def fitness_func(solution, solution_idx):
    global generation
    global step
    global population_avg_list
    global population_min_avg_list
    global even_list
    global uneven_list

    population_avg_list = []
    population_min_avg_list = []

    os.chdir(os.environ['WORKING_DIR']+'/tmp')

    dir_name = str(uuid.uuid4())
    os.mkdir(os.environ['WORKING_DIR']+'/tmp/'+dir_name)
    os.chdir(os.environ['WORKING_DIR']+'/tmp/'+dir_name)

    os.environ['TMP_DIR']=os.getcwd()

    os.system("cp -r ../../protein_* .")

    for n_prot in dir_list:

        os.chdir(os.environ['TMP_DIR']+f'/protein_{n_prot}')
        os.environ['OUTPUT_DIR']=os.environ['TMP_DIR']+f'/protein_{n_prot}/output'

        # Obtain ligand and protein names
        for files in glob.glob("*.xyz"):
            file_list = files.split('_c.')
            name_ligand = file_list[0]

        for files in glob.glob("*.pdb"):
            file_list = files.split('.')
            name_protein = file_list[0]

        if iv.var.scale_factor != None:
            npts = box_size_func('*.sdf', 0.375, iv.var.scale_factor)
        
        if iv.var.box_size != None: 
            npts = iv.var.box_size
        
        if iv.var.step_wise == True:    
            try:
                npts = box_size_func('*.sdf', 0.375, box_size[step])
            except IndexError:
                npts = box_size_func('*.sdf', 0.375, box_size[-1])

        ##### AutoDock ##### 
        os.chdir(os.environ['OUTPUT_DIR'])

        if os.path.isdir('docking') == False:
            os.mkdir('docking')
            os.chdir('docking')
        else:
            os.chdir('docking')

        e = open('energy','r')
        lines = [line.split() for line in e]
        energy = lines[0][4]

        dock = docking_centre('coordinates')

        docking_func(solution, name_ligand, name_protein, energy, dock, npts, solution_idx)

        ##### Fitness function ######
        avg_list, min_list, output = rmsd_func(name_ligand, n_prot)

        population_avg_list.append(avg_list)
        population_min_avg_list.append(min_list)


    '''Compare previous generations with each. If the average of the solutions between each two generations is below a certain threshold
    then the size of the box will be increased.
    ''' 
    if iv.var.step_wise == True:  
        if generation % 2 == 0:
            "Even generations"
            even_list = np.append(even_list,solution)
            
            if generation != 0:
                difference = [np.abs(i[0]-i[1]) for i in zip(even_list,uneven_list)]
                
                if np.mean(np.array(difference)) < 0.5:
                    step+=1
                    #print('Item of boxsize list is now item {}'.format(step))

                even_list = np.zeros([1,12])
                    

        if generation % 2 != 0:
            "Uneven generations"
            uneven_list = np.append(uneven_list,solution)

            difference = [np.abs(i[0]-i[1]) for i in zip(even_list,uneven_list)]
            
            if np.mean(np.array(difference)) < 0.5:
                step+=1
                #print('Item of boxsize list is now item {}'.format(step))
            
            uneven_list = np.zeros([1,12])

    os.chdir(os.environ['WORKING_DIR'])
    """ Method to Calculate Fitness

    One protein:
                avg_output: takes the average of the RMSD of the poses
                minimum_rmsd: takes the minimum of the RMSD of the poses

    Multiple protein:
                generation_avg: takes the average of the RMSD of the poses of each protein and averages once more
                generation_min_avg: takes the minimum of the RMSD of the poses of each protein and averages
    """
    if len(dir_list) > 1:
        sum_population = sum(flatten(population_avg_list))
        sum_min_population = sum(flatten(population_min_avg_list))

        population_avg = sum_population/ len(population_avg_list)
        population_min_avg = sum_min_population / len(population_min_avg_list)

        fitness = 1 / np.abs(population_avg - desired_output)            

        output = [f"-------------------------------------------     PARENT  {solution_idx}      ---------------------------------------------\n"]
        output.append(f"Average RMSD of parent {solution_idx} in generation {generation}: {population_avg:.4}\n")
        output.append(f"Average RMSD of the lowest RMSDs of parent {solution_idx} in generation {generation}: {population_min_avg:.4}\n")
        output.append("-----------------------------------------------------------------------------------------------------------\n")

        with open('parameter_history', 'a') as f:
            f.write('Parameters generation {}         :  '.format(str(generation))+'  '.join(format(solution[x], ">10.5f") for x in range(0,len(solution)))+'| {:>10.5f}  {:>10.5f}     {:>10.5f}\n'.format(fitness,population_avg, population_min_avg))
    
    else:
        avg_output = np.mean(np.array((avg_list)))
        min_avg_rmsd = np.mean(np.array(min_list))

        fitness = 1 / np.abs(avg_output - desired_output)

        output = [f"-------------------------------------------     PARENT  {solution_idx}      ---------------------------------------------\n"]
        output.append(f"Average RMSD of parent {solution_idx} in generation {generation}: {avg_output:.4}\n")
        output.append(f"Average RMSD of the lowest RMSDs of parent {solution_idx} in generation {generation}: {min_avg_rmsd:.4}\n")
        output.append("-----------------------------------------------------------------------------------------------------------\n")

        with open('parameter_history', 'a') as f:
            f.write('Parameters generation {}         :  '.format(str(generation))+'  '.join(format(solution[x], ">10.5f") for x in range(0,len(solution)))+'| {:>10.5f}  {:>10.5f}     {:>10.5f}\n'.format(fitness, avg_output, min_avg_rmsd))

    generation+=1
    shutil.rmtree(os.environ['TMP_DIR'],ignore_errors=True)

    print(''.join(output))

    return fitness
    

def fitness_wrapper(solution, solution_idx):
    return fitness_func(solution, solution_idx)

def plot_each_protein(protein_numbers):
    mpl.rcParams['figure.dpi'] = 400

    plt.rcParams.update({'font.size': 12})
    plt.rcParams['font.family'] = 'DeJavu Serif'
    plt.rcParams['font.serif'] = ['Palatino']

    first_gen = []
    last_gen = []

    protein_number = protein_numbers

    first_gen_dict = {}
    last_gen_dict = {}

    first_protein_rmsd_dict = {}
    last_protein_rmsd_dict = {}

    for x in range(0,len(protein_number)):
        first_gen_dict[x] = protein_number[x]
        first_protein_rmsd_dict[protein_number[x]] = []
        
        first_gen_dict[x] = os.environ['WORKING_DIR']+'/first_gen/protein_'+str(first_gen_dict[x])
        
        protein = open(first_gen_dict[x])
        protein_rmsd = [line.split() for line in protein]
        
        for i in range(0,len(protein_rmsd)):
            first_gen.append(float(protein_rmsd[i][5]))
            first_protein_rmsd_dict[protein_number[x]].append(float(protein_rmsd[i][5]))
            
        last_gen_dict[x] = protein_number[x]
        last_protein_rmsd_dict[protein_number[x]] = []
        
        last_gen_dict[x] = os.environ['WORKING_DIR']+'/last_gen/protein_'+str(last_gen_dict[x])
        
        protein = open(last_gen_dict[x])
        protein_rmsd = [line.split() for line in protein]
        
        for i in range(0,len(protein_rmsd)):
            last_gen.append(float(protein_rmsd[i][5]))
            last_protein_rmsd_dict[protein_number[x]].append(float(protein_rmsd[i][5]))

        mu_first = statistics.mean(first_protein_rmsd_dict[protein_number[x]])
        sigma_first = statistics.stdev(first_protein_rmsd_dict[protein_number[x]])
        
        mu_last = statistics.mean(last_protein_rmsd_dict[protein_number[x]])
        sigma_last = statistics.stdev(last_protein_rmsd_dict[protein_number[x]])
        
        label_first = [f'First Generation: $\mu$ = {mu_first:1.2f} $\sigma$ = {sigma_first:1.2f}']
        label_last = [f'Last Generation: $\mu$ = {mu_last:1.2f} $\sigma$ = {sigma_last:1.2f}']

        bins = np.linspace(0, 15, 100)

        fig, ax = plt.subplots(1, figsize=(6,4), sharex=True)

        ax.hist(first_protein_rmsd_dict[protein_number[x]], density=False, color='b', alpha=0.5, bins=bins, label=label_first)# density=False would make counts
        ax.hist(last_protein_rmsd_dict[protein_number[x]], density=False, color='r', alpha=0.5, bins =bins, label=label_last)

        ax.set_title('Protein {}'.format(protein_number[x]))

        ax.set_xticks(range(0,15,2))
        
        ax.xaxis.set_major_formatter(FormatStrFormatter('%i'))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%i'))

        ax.set_xlim(0,15)
        ax.set_ylabel('Conformations (N)')
        ax.set_xlabel('RMSD (Å)')

        plt.legend(loc='upper right')
    
        plt.tight_layout()
        plt.savefig(f"protein_"+str(protein_number[x])+".png", bbox_inches='tight')
        plt.close()

def plot_total_conformations():
    mpl.rcParams['figure.dpi'] = 400

    plt.rcParams.update({'font.size': 12})
    plt.rcParams['font.family'] = 'DeJavu Serif'
    plt.rcParams['font.serif'] = ['Palatino']

    first_gen_all = []
    last_gen_all = []

    first_gen = open(os.environ['WORKING_DIR']+'/first_gen/all_conf_first_gen')
    first_rmsd = [line.split() for line in first_gen]
        
    for i in range(0,len(first_rmsd)):
        first_gen_all.append(float(first_rmsd[i][5]))

    last_gen = open(os.environ['WORKING_DIR']+'/last_gen/all_conf_last_gen')
    last_rmsd = [line.split() for line in last_gen]
        
    for i in range(0,len(last_rmsd)):
        last_gen_all.append(float(last_rmsd[i][5]))

    mu_first = statistics.mean(first_gen_all)
    sigma_first = statistics.stdev(first_gen_all)
    
    mu_last = statistics.mean(last_gen_all)
    sigma_last = statistics.stdev(last_gen_all)
    
    label_first = [f'First Generation: $\mu$ = {mu_first:1.2f} $\sigma$ = {sigma_first:1.2f}']
    label_last = [f'Last Generation: $\mu$ = {mu_last:1.2f} $\sigma$ = {sigma_last:1.2f}']

    bins = np.linspace(0, 15, 100)

    fig, ax = plt.subplots(1, figsize=(6,4), sharex=True)
    
    ax.hist(first_gen_all, density=False, color='b', alpha=0.5, bins=bins, label=label_first)# density=False would make counts
    ax.hist(last_gen_all, density=False, color='r', alpha=0.5, bins =bins, label=label_last)
    
    ax.legend(loc="upper right")
    ax.set_title('All Conformations')
    
    ax.set_xlabel('RMSD (Å)')
    ax.set_ylabel('Conformations(N)')
    
    ax.set_xticks(range(0,15,2))

    ax.xaxis.set_major_formatter(FormatStrFormatter('%i'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%i'))

    plt.legend(loc='upper right')

    plt.tight_layout()
    plt.savefig(f"all_conformations.png", bbox_inches='tight')
    plt.close()

def plot_parameters(metal):
    mpl.rcParams['figure.dpi'] = 400

    plt.rcParams.update({'font.size': 12})
    plt.rcParams['font.family'] = 'DeJavu Serif'
    plt.rcParams['font.serif'] = ['Palatino']

    data_set = open(os.environ['WORKING_DIR']+'/parameter_history','r')
    data_set_line = [ lines.split() for lines in data_set]

    fitness = []
    rmsd_min = []

    r_OA = []
    e_OA = []

    r_SA = []
    e_SA = []

    r_HD = []
    e_HD = []

    r_NA = []
    e_NA = []

    r_N = []
    e_N = []

    r_M = []
    e_M = []
    

    for i in range(1,len(data_set_line)):
        r_OA.append(float(data_set_line[i][4]))
        e_OA.append(float(data_set_line[i][5]))
        
        r_SA.append(float(data_set_line[i][6]))
        e_SA.append(float(data_set_line[i][7]))
        
        r_HD.append(float(data_set_line[i][8]))
        e_HD.append(float(data_set_line[i][9]))
        
        r_NA.append(float(data_set_line[i][8]))
        e_NA.append(float(data_set_line[i][9]))
        
        r_N.append(float(data_set_line[i][10]))
        e_N.append(float(data_set_line[i][11]))
        
        r_M.append(float(data_set_line[i][12]))
        e_M.append(float(data_set_line[i][13]))
        
        fitness.append(float(data_set_line[i][16]))
        rmsd_min.append(float(data_set_line[i][18]))

    x_axis = range(0,len(r_OA))

    fig, ax = plt.subplots(6,2, figsize=(15,15), sharex=True)

    ax[0][0].plot(x_axis,r_OA, c='b')
    ax[0][1].plot(x_axis,e_OA, c='b')

    ax[1][0].plot(x_axis,r_SA, c='b')
    ax[1][1].plot(x_axis,e_SA, c='b')

    ax[2][0].plot(x_axis,r_HD, c='b')
    ax[2][1].plot(x_axis,e_HD, c='b')

    ax[2][0].plot(x_axis,r_NA, c='b')
    ax[2][1].plot(x_axis,e_NA, c='b')

    ax[3][0].plot(x_axis,r_N, c='b')
    ax[3][1].plot(x_axis,e_N, c='b')

    ax[4][0].plot(x_axis,r_N, c='b')
    ax[4][1].plot(x_axis,e_N, c='b')

    ax[5][0].plot(x_axis,r_M, c='b')
    ax[5][1].plot(x_axis,e_M, c='b')

    #########################################################
    ax[0][0].set_ylabel('r {}-OA (Å)'.format(metal))
    ax[0][1].set_ylabel('$\epsilon$ {}-OA (kcal/mol)'.format(metal))

    ax[1][0].set_ylabel('r {}-SA (Å)'.format(metal))
    ax[1][1].set_ylabel('$\epsilon$ {}-SA (kcal/mol)'.format(metal))

    ax[2][0].set_ylabel('r {}-HD (Å)'.format(metal))
    ax[2][1].set_ylabel('$\epsilon$ {}-HD (kcal/mol)'.format(metal))

    ax[3][0].set_ylabel('r {}-NA (Å)'.format(metal))
    ax[3][1].set_ylabel('$\epsilon$ {}-NA (kcal/mol)'.format(metal))

    ax[4][0].set_ylabel('r {}-N (Å)'.format(metal))
    ax[4][1].set_ylabel('$\epsilon$ {}-N (kcal/mol)'.format(metal))
        
    ax[5][0].set_ylabel('r {}-{} (Å)'.format(metal,metal))
    ax[5][1].set_ylabel('$\epsilon$ {}-{} (kcal/mol)'.format(metal,metal))

    ax[5][0].set_xlabel('N parents')
    ax[5][1].set_xlabel('N parents')

    # for i in range(0,6):
    #     for j in range(0,2):
    #         ax[i][j].set_xlim(0,500)

    plt.tight_layout()
    plt.savefig(f"parameters.png", bbox_inches='tight')
    plt.close()

    fig, ax = plt.subplots(1,2, figsize=(15,5), sharex=True)

    ax[0].plot(x_axis,fitness, c='b', label='Data set')

    ax[1].plot(x_axis,rmsd_min, c='b', label='{}'.format(metal))

    ax[0].set_ylabel('Fitness Function')
    ax[1].set_ylabel('RMSD (Å)')

    ax[0].set_xlabel('N parentes')
    ax[1].set_xlabel('N parents')

    plt.legend()
    plt.tight_layout()
    plt.savefig(f"RMSD_fitness.png", bbox_inches='tight')
    plt.close()


class PooledGA(pygad.GA):

    def cal_pop_fitness(self):
        global pool

        pop_fitness = pool.starmap(fitness_wrapper, zip(self.population, list(range(0,len(self.population)))))
        pop_fitness = np.array(pop_fitness)
        return pop_fitness

def train_GA():

    iv.insert_arguments()

    os.environ['WORKING_DIR']=os.getcwd()

    if os.path.exists('parameter_history'):
        os.remove('parameter_history')

    if os.path.isdir(os.environ['WORKING_DIR']+'/first_gen'):
        shutil.rmtree(os.environ['WORKING_DIR']+'/first_gen', ignore_errors=True)
        os.mkdir(os.environ['WORKING_DIR']+'/first_gen')
    else:
        os.mkdir(os.environ['WORKING_DIR']+'/first_gen')

    if os.path.isdir(os.environ['WORKING_DIR']+'/last_gen'):
        shutil.rmtree(os.environ['WORKING_DIR']+'/last_gen', ignore_errors=True)
        os.mkdir(os.environ['WORKING_DIR']+'/last_gen')
    else:
        os.mkdir(os.environ['WORKING_DIR']+'/last_gen')

    # Global variables
    generation = 0
    step = 0

    with open('parameter_history', 'a') as f:
        f.write("All old solutions are           :     r_OA        e_OA        r_SA        e_SA        r_HD        e_HD        r_NA        e_NA        r_N         e_N         r_"+iv.var.metal_symbol+"_"+iv.var.metal_symbol+"     e_"+iv.var.metal_symbol+"_"+iv.var.metal_symbol+"|    fitness    RMSD_AVG   RMSD_MIN_AVG\n")

    # Make list of the protein numbers to iterate over
    dir_list = os.listdir(os.getcwd())
    dir_list = [str(i).replace('protein_','') for i in dir_list]
    dir_list = [int(i) for i in dir_list if convertible(i)]
    dir_list = sorted(dir_list)

    # # Parameters 
    desired_output = 0

    fitness_function = fitness_func

    # Scale factor list
    box_size = [1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0]

    num_generations = iv.var.num_generations
    num_parents_mating = iv.var.num_parents_mating

    # initial_population = None

    sol_per_pop = iv.var.sol_per_pop
    num_genes = 12
    gene_space=[{'low': 1, 'high': 3},{'low': 0, 'high': 25},
                {'low': 1, 'high': 3},{'low': 0, 'high': 25},
                {'low': 1, 'high': 3},{'low': 0, 'high': 25},
                {'low': 1, 'high': 3},{'low': 0, 'high': 25},
                {'low': 1, 'high': 3},{'low': 0, 'high': 25},
                {'low': 1, 'high': 3},{'low': 0, 'high': 25}]

    parent_selection_type = iv.var.parent_selection_type
    keep_parents = iv.var.keep_parents
    k_tournament = iv.var.k_tournament

    crossover_type = iv.var.crossover_type
    crossover_probability = iv.var.crossover_prob

    mutation_probability = iv.var.mutation_prob
    mutation_type = mutation_func

    # Create Class
    ga_instance = PooledGA(num_generations=num_generations,
                           num_parents_mating=num_parents_mating,
                           fitness_func=fitness_function,
                           sol_per_pop=sol_per_pop,
                           num_genes=num_genes,
                           gene_space=gene_space,
                           parent_selection_type=parent_selection_type,
                           keep_parents=keep_parents,
                           K_tournament=k_tournament,
                           crossover_type=crossover_type,
                           crossover_probability=crossover_probability,
                           mutation_type=mutation_func,
                           save_solutions=True)

    os.mkdir(os.environ['WORKING_DIR']+'/tmp')
    os.chdir(os.environ['WORKING_DIR']+'/tmp')

    os.environ['TMP_DIR']=os.getcwd()

    even_list = np.zeros([1,12])
    uneven_list = np.zeros([1,12])

    with Pool(processes=sol_per_pop) as pool:
        ga_instance.run()

        solution, solution_fitness, solution_idx = ga_instance.best_solution()
        print("Parameters of the best solution : {solution}\n".format(solution=solution))
        print("Fitness value of the best solution = {solution_fitness}\n".format(solution_fitness=solution_fitness))

    shutil.rmtree(os.environ['WORKING_DIR']+'/tmp',ignore_errors=True)


    if os.path.isdir(os.environ['WORKING_DIR']+'/figures'):
        shutil.rmtree(os.environ['WORKING_DIR']+'/figures', ignore_errors=True)
        os.mkdir(os.environ['WORKING_DIR']+'/figures')
    else:
        os.mkdir(os.environ['WORKING_DIR']+'/figures')

    os.chdir(os.environ['WORKING_DIR']+'/figures')

    plot_each_protein(dir_list)
    plot_total_conformations()
    plot_parameters(iv.var.metal_symbol)
