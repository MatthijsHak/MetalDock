#!/usr/bin/env python3 

import argparse
import glob,sys,os
import subprocess
import numpy as np
import pygad
import uuid
import shutil
import random

from multiprocessing import Pool
from distutils.dir_util import copy_tree


from openbabel import pybel as py
from openbabel import openbabel as ob

import docking as dock
import input_variables as iv
import environment_variables

ob.obErrorLog.SetOutputLevel(0)

def convertible(v):
    try:
        int(v)
        return True
    except (TypeError, ValueError):
        return False

def mutation_func(offspring, ga_instance):
    '''
    Mutation function that adds 10 percent or substracts 10 percent of the parameter value
    '''
    print(mutation_probability)
    for chromosome_idx in range(0,len(offspring[0])):
        random_number_1 = random.uniform(0,1)
        random_number_2 = random.uniform(0,1)

        if random_number_1 < mutation_probability:
            if random_number_2 <= 0.5:
                mutated_gene = offspring[0][chromosome_idx]*1.1
            else:
                mutated_gene = offspring[0][chromosome_idx]*0.9

            offspring[0][chromosome_idx] = mutated_gene

    return offspring

def fitness_func(solution, solution_idx):
    global parent
    global generation

    generation_avg_list = []
    generation_min_avg_list = []

    os.chdir(os.environ['WORKING_DIR']+'/tmp')

    dir_name = str(uuid.uuid4())
    os.mkdir(os.environ['WORKING_DIR']+'/tmp/'+dir_name)
    os.chdir(os.environ['WORKING_DIR']+'/tmp/'+dir_name)

    os.environ['TMP_DIR']=os.getcwd()
    os.system("cp -r ../../protein* .")

    for n_prot in dir_list:

        os.chdir(os.environ['WORKING_DIR']+f'/protein_{n_prot}')

        for files in glob.glob("*.xyz"):
            file_list = files.split('.')
            name_ligand = file_list[0]

        for files in glob.glob("*.pdb"):
            file_list = files.split('.')
            name_protein = file_list[0]

        os.environ['OUTPUT_DIR']=os.environ['TMP_DIR']+f'/protein_{n_prot}/output'

        ##### AutoDock ##### 
        os.chdir(os.environ['OUTPUT_DIR']+'/file_prep')

        dock_site = open('coordinates','r')
        coord = [line.split() for line in dock_site]

        dock_x = coord[0][0]
        dock_y = coord[0][1]
        dock_z = coord[0][2]

        os.chdir(os.environ['OUTPUT_DIR'])

        if os.path.isdir('docking') == False:
            os.mkdir('docking')
            os.chdir('docking')
        else:
            os.chdir('docking')

        dock.randomize_translation_rotation(name_ligand+'.pdbqt')
        os.system('mv docking.pdbqt '+name_ligand+'.pdbqt')

        os.system('cp '+os.environ['WORKING_DIR']+'/'+iv.var.parameter_file+' .')

        # insert solution for R and epsilon for H-bond
        os.system(r'''awk '{ if ($2 == "RU" || $2 == "Ru") ($7 = '''+str(solution[10])+''') && ($8 = '''+str(solution[11])+'''); print $0}' '''+iv.var.parameter_file+''' > file_1''')
        os.system(r'''awk '{ if ($2 == "RU" || $2 == "Ru") printf "%-8s %-3s %7s %8s %8s %9s %4s %4s %2s %3s %3s %2s\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12; else print $0}' file_1 > '''+iv.var.parameter_file)

        #add_to_dat_file():
        dat = open(''+iv.var.parameter_file+'', 'a')
        dat.write('nbp_r_eps '+str(solution[0])+'   '+str(solution[1])+'    12 6 OA '+iv.var.metal_symbol+'\n')
        dat.write('nbp_r_eps '+str(solution[2])+'   '+str(solution[3])+'    12 6 SA '+iv.var.metal_symbol+'\n')
        dat.write('nbp_r_eps '+str(solution[4])+'   '+str(solution[5])+'    12 6 HD '+iv.var.metal_symbol+'\n')
        dat.write('nbp_r_eps '+str(solution[6])+'   '+str(solution[7])+'    12 6 NA '+iv.var.metal_symbol+'\n')
        dat.write('nbp_r_eps '+str(solution[8])+'   '+str(solution[9])+'    12 6  N '+iv.var.metal_symbol+'\n')
        dat.close()


        #create_gpf():
        gpf = open('clean_'+name_protein+'.gpf', 'a')
        gpf.write('nbp_r_eps '+str(solution[0])+'   '+str(solution[1])+'    12 6 OA '+iv.var.metal_symbol+'\n')
        gpf.write('nbp_r_eps '+str(solution[2])+'   '+str(solution[3])+'    12 6 SA '+iv.var.metal_symbol+'\n')
        gpf.write('nbp_r_eps '+str(solution[4])+'   '+str(solution[5])+'    12 6 HD '+iv.var.metal_symbol+'\n')
        gpf.write('nbp_r_eps '+str(solution[6])+'   '+str(solution[7])+'    12 6 NA '+iv.var.metal_symbol+'\n')
        gpf.write('nbp_r_eps '+str(solution[8])+'   '+str(solution[9])+'    12 6  N '+iv.var.metal_symbol+'\n')
        gpf.close()

        #autogrid()
        os.system(os.environ['AUTODOCK']+'/autogrid4 -p clean_'+name_protein+'.gpf > /dev/null 2>&1')

        #create_dpf()
        os.system(os.environ['PYTHON_2']+" "+os.environ['MGLTOOLS']+"/prepare_dpf42.py -l "+name_ligand+".pdbqt -r clean_"+name_protein+".pdb -p parameter_file="+iv.var.parameter_file+" > /dev/null 2>&1")

        #autodock()
        os.system(os.environ['AUTODOCK']+'/autodock4 -p '+name_ligand+'_clean_'+name_protein+'.dpf > /dev/null 2>&1')

        #write_all_conformations()
        os.system(os.environ['PYTHON_2']+" "+os.environ['MGLTOOLS']+"/write_conformations_from_dlg.py -d "+name_ligand+"_clean_"+name_protein+".dlg")

        #pdb = next(py.readfile('pdb','ref.pdb'))
        #pdb.write('xyz','ref.xyz',overwrite=True)

        rmsd_avg = []
        rmsd_list = []

        os.system(os.environ['OBABEL']+" -isdf output.sdf -oxyz normalize.xyz -d > normalize.xyz")

        rmsd_normalize = float(subprocess.getoutput([os.environ['PYTHON_3']+' '+os.environ['LIB_DIR']+'/calculate_rmsd.py ref.xyz normalize.xyz --reorder']))
        rmsd_list.append("RMSD between reference ligand and quantum optimized structure: %.4f" % rmsd_normalize)

        i = 1
        while os.path.exists(name_ligand+"_%i.pdbqt" % i):
            os.system(os.environ['OBABEL']+" -ipdbqt "+name_ligand+"_{}.pdbqt".format(i)+" -oxyz "+name_ligand+"_{}.xyz".format(i)+" -d > "+name_ligand+"_{}.xyz".format(i))

            rmsd_non_rotate = float(subprocess.getoutput([os.environ['PYTHON_3']+' '+os.environ['LIB_DIR']+'/calculate_rmsd.py ref.xyz '+name_ligand+'_{}.xyz'.format(i)+' --reorder --rotation none --translation none']))
            rmsd = rmsd_non_rotate

            rmsd_list.append("RMSD for Conformation %i = %.4f"% (i, rmsd))
            rmsd_avg.append(rmsd)
            i += 1

        for j in range(0,len(rmsd_list)):
            print(rmsd_list[j])

        avg_output = (sum(rmsd_avg) / len(rmsd_avg))
        generation_avg_list.append(avg_output)
        print(f"The average RMSD of protein {n_prot}: {avg_output:.4}")

        minimum_rmsd = min(rmsd_avg)
        generation_min_avg_list.append(minimum_rmsd)
        print(f"The lowest RMSD value of the poses of protein {n_prot}: {minimum_rmsd:.4}")
        print("-----------------------------------------------------------------------------------------------------------")


    os.chdir(os.environ['WORKING_DIR'])
    """Method to Calculate Fitness

    One protein:
                avg_output: takes the average of the RMSD of the poses
                minimum_rmsd: takes the minimum of the RMSD of the poses

    Multiple protein:
                generation_avg: takes the average of the RMSD of the poses of each protein and averages once more
                generation_min_avg: takes the minimum of the RMSD of the poses of each protein and averages
    """
    if len(dir_list) > 1:
            generation_avg = (sum(generation_avg_list) / len(generation_avg_list))
            generation_min_avg= (sum(generation_min_avg_list) / len(generation_min_avg_list))

    method = generation_avg

    fitness = 1 / np.abs(method - desired_output)


    if len(dir_list) > 1:
        print(f"The average RMSD of all proteins in this parameter set in generation {generation}: {generation_avg:.4}")
        print(f"The average RMSD of the lowest conformation of each ligand in this parameter set in generation {generation}: {generation_min_avg:.4}")
        print("-----------------------------------------------------------------------------------------------------------")

    with open('parameter_history', 'a') as f:
        f.write('Parameters generation {}         :  '.format(str(generation))+'  '.join(format(solution[x], ">10.5f") for x in range(0,len(solution)))+'| {:>10.5f}  {:>10.5f}     {:>10.5f}\n'.format(fitness,generation_avg, generation_min_avg))

    generation+=1
    shutil.rmtree(os.environ['TMP_DIR'],ignore_errors=True)

    return fitness


def fitness_wrapper(solution):
    return fitness_func(solution, 0)


class PooledGA(pygad.GA):

    def cal_pop_fitness(self):
        global pool

        pop_fitness = pool.map(fitness_wrapper, self.population)
        pop_fitness = np.array(pop_fitness)
        return pop_fitness

if __name__=='__main__':

    iv.insert_arguments()

    os.environ['WORKING_DIR']=os.getcwd()

    if os.path.exists('parameter_history'):
        os.remove('parameter_history')

    parent = 1
    generation = 0

    with open('parameter_history', 'a') as f:
        f.write("All old solutions are           :     r_OA        e_OA        r_SA        e_SA        r_HD        e_HD        r_NA        e_NA        r_N         e_N         r_Ru_Ru     e_Ru_Ru|    fitness    RMSD_AVG   RMSD_MIN_AVG\n")

    # Make list of the protein numbers to iterate over
    dir_list = os.listdir(os.getcwd())
    dir_list = [str(i).replace('protein_','') for i in dir_list]
    dir_list = [int(i) for i in dir_list if convertible(i)]
    dir_list = sorted(dir_list)

    # Parameters 
    desired_output = 0

    fitness_function = fitness_func

    num_generations = iv.var.num_generations
    num_parents_mating = iv.var.num_parents_mating

    initial_population = None

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
    K_tournament = iv.var.K_tournament

    crossover_type = iv.var.crossover_type
    crossover_probability = iv.var.crossover_prob

    #mutation_type = iv.var.mutation_type
    mutation_probability = iv.var.mutation_prob
    mutation_percent_genes = iv.var.mutation_percent

    # Create Class
    #if mutation_probability != None:
    ga_instance = PooledGA(num_generations=num_generations,
                           num_parents_mating=num_parents_mating,
                           fitness_func=fitness_function,
                           initial_population=initial_population,
                           sol_per_pop=sol_per_pop,
                           num_genes=num_genes,
                           gene_space=gene_space,
                           parent_selection_type=parent_selection_type,
                           keep_parents=keep_parents,
                           K_tournament=K_tournament,
                           crossover_type=crossover_type,
                           crossover_probability=crossover_probability,
                           mutation_type=mutation_func,
                           save_solutions=True)

    os.mkdir(os.environ['WORKING_DIR']+'/tmp')

    with Pool(processes=sol_per_pop) as pool:
        ga_instance.run()

    solution, solution_fitness, solution_idx = ga_instance.best_solution()
    print("Parameters of the best solution : {solution}".format(solution=solution))
    print("Fitness value of the best solution = {solution_fitness}".format(solution_fitness=solution_fitness))

    shutil.rmtree(os.environ['WORKING_DIR']+'/tmp',ignore_errors=True)
