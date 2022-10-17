#!/usr/bin/env python3 

import os,sys,glob
import subprocess
import pygad
import uuid
import shutil
import random
import environment_variables

import numpy as np
import multiprocessing as mp 
import prepare_dock as d
import figures as fig
import optim as opt

from scipy.stats import rankdata
from multiprocessing import Pool
from parser import Parser

input_dir = os.getcwd()
output_dir = f'{input_dir}/output'

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
    mutation_step = 100 
    mutation_positive = [x * 0.1 for x in range(1, mutation_step)]
    mutation_negative = [x * 0.001 for x in range(1, mutation_step)]

    #Adaptive Mutation
    mutation_probability_list = [ x  for x in np.linspace(0, 1, num=par.sol_per_pop)]

    rank_list = rankdata(ga_instance.last_generation_fitness)

    offspring_length = par.sol_per_pop - par.keep_par

    for i in range(0,offspring_length):
        mutation_probability = mutation_probability_list[int(rank_list[i])-1]

        for chromosome_idx in range(0,len(ga_instance.gene_space)):
            random_number_1 = random.uniform(0,1)
            random_number_2 = random.uniform(0,1)
            random_int = random.randrange(0, len(mutation_positive))

            if random_number_1 < mutation_probability:
                if random_number_2 <= 0.5:
                    mutated_gene = offspring[i][chromosome_idx] * mutation_positive[random_int]
                else:
                    mutated_gene = offspring[i][chromosome_idx] * mutation_negative[random_int]

                if ga_instance.gene_space[chromosome_idx]['low'] <= mutated_gene <= ga_instance.gene_space[chromosome_idx]['high']:
                    pass
                else:
                    mutated_gene = offspring[i][chromosome_idx]
            
            else:
                offspring[i] = offspring[i]

    return offspring

# def docking_centre(coordinate_file):
#     dock_site = open(coordinate_file,'r')
#     coord = [line.split() for line in dock_site]

#     if is_float(coord[0][0]) == True:
#         dock = [coord[0][0], coord[0][1], coord[0][2]]
#     else:
#         dock = [coord[0][1], coord[0][2], coord[0][3]]

#     return dock

# def rmsd_func(name_ligand, n_prot):
#     rmsd_avg = []
#     rmsd_list = []
#     avg_list = []
#     min_list = []
#     rmsd_print_list = []
                
#     output = [f"-------------------------------------------     PROTEIN {n_prot}      --------------------------------------------\n"]
              
#     i = 1
#     while os.path.exists(name_ligand+"_%i.pdbqt" % i):
#         subprocess.call(os.environ['OBABEL']+" -ipdbqt "+name_ligand+"_{}.pdbqt".format(i)+" -oxyz "+name_ligand+"_{}.xyz".format(i)+" -d > "+name_ligand+"_{}.xyz".format(i), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) 

#         rmsd_non_rotate = float(subprocess.getoutput([os.environ['PYTHON_3']+' '+os.environ['DOCK_LIB_DIR']+'/calculate_rmsd.py ref.xyz '+name_ligand+'_{}.xyz'.format(i)+' --reorder --rotation none --translation none']))
#         rmsd = rmsd_non_rotate

#         rmsd_print_list.append("RMSD for Conformation %i = %.4f\n"% (i, rmsd))
#         rmsd_list.append(rmsd)
#         i += 1

#     for j in range(0,len(rmsd_print_list)):
#         output.append(rmsd_print_list[j])

#         if gen == 0:
#             first_gen = open(f'{output_dir}/first_gen/all_conf_first_gen','a')
#             first_gen.write(rmsd_print_list[j])

#             protein = open(f'{output_dir}/first_gen/protein_{n_prot}','a')
#             protein.write(rmsd_print_list[j])

#         if gen == par.ga_num_generations+1:   
#             last_gen = open(f'{output_dir}/last_gen/all_conf_last_gen','a')
#             last_gen.write(rmsd_print_list[j])
            
#             protein = open(f'{output_dir}/last_gen/protein_{n_prot}','a')
#             protein.write(rmsd_print_list[j])

#     avg_output = np.mean(rmsd_list)
#     avg_list.append(avg_output)
#     output.append(f"Average RMSD: {avg_output:.4}\n")

#     minimum_rmsd = min(rmsd_list)
#     min_list.append(minimum_rmsd)
#     output.append(f"Lowest RMSD: {minimum_rmsd:.4}\n")

#     stdv_rmsd = np.std(rmsd_list)
#     output.append(f"Standard Deviation RMSD: {stdv_rmsd:.4}\n")

#     var_rmsd = np.var(rmsd_list)
#     output.append(f"Variance RMSD: {var_rmsd:.4}\n")
#     output.append("-----------------------------------------------------------------------------------------------------------\n")

#     return avg_list, min_list, print(''.join(output))
    


def fitness_func(solution, solution_idx):
    global step
    global gen
    global population_avg_list
    global population_min_avg_list
    global even_list
    global uneven_list

    # Parameters 
    desired_output = 0

    population_avg_list = []
    population_min_avg_list = []

    os.chdir(f'{tmp_dir}')

    dir_name = str(uuid.uuid4())
    os.mkdir(f'{tmp_dir}/'+dir_name)
    os.chdir(f'{tmp_dir}/'+dir_name)

    tmp_dir_GA=os.getcwd()

    # os.system("cp -r ../../protein_* .")

    for n_prot in dir_list:
        os.chdir(f'{tmp_dir}/'+dir_name)
        os.mkdir(f'protein_{n_prot}')
        os.chdir(f'{tmp_dir}/'+dir_name+f'/protein_{n_prot}')

        os.system(f"cp -r  ../../../../data_set/protein_{n_prot}/output/docking .")
        os.chdir(f'{tmp_dir}/'+dir_name+f'/protein_{n_prot}/docking')
        output_temp_dir=f'{tmp_dir}/'+dir_name+f'/protein_{n_prot}/docking'

        # Obtain ligand and protein names
        for files in glob.glob("*_1.pdbqt"):
            file_list = files.split('_1.pdbqt')
            name_ligand = file_list[0]

        for files in glob.glob("clean_*.pdb"):
            file_list = files.split('clean_')
            file_list = file_list[1].split('.')
            name_protein = file_list[0]


        if par.scale_factor > 0:
            npts = d.box_size_func('ref.xyz', par.metal_symbol, 0.375, par.scale_factor)
        
        if par.box_size > 0: 
            npts = [par.box_size, par.box_size, par.box_size]
        
        if par.step_wise == True:    
            try:
                npts = d.box_size_func('ref.xyz', par.metal_symbol, 0.375, box_size[step])
            except IndexError:
                npts = d.box_size_func('ref.xyz', par.metal_symbol, 0.375, box_size[-1])

        ##### AutoDock ##### 
        e = open(f'{output_temp_dir}/energy','r')
        lines = [line.split() for line in e]
        energy = lines[0][4]

        dock = d.get_coordinates(f'{output_temp_dir}/ref.xyz',par.metal_symbol)

        subprocess.call([f'cp {input_dir}/'+par.parameter_file+' .'], shell=True)
        d.docking_func(solution, par.parameter_file, par.metal_symbol, name_ligand, name_protein, energy, dock, npts, par.num_poses, par.dock_algorithm, par.random_pos, par.ga_dock, par.sa_dock)

        ##### Fitness function ######
        avg_list, min_list, output = d.rmsd_func(name_ligand, n_prot, gen, output_dir, num_gen=par.ga_num_generations, train=True)

        population_avg_list.append(avg_list)
        population_min_avg_list.append(min_list)

        subprocess.call([f'rm -r {tmp_dir}/'+dir_name+f'/protein_{n_prot}'], shell=True)


    '''Compare previous generations with each. If the average of the solutions between each two generations is below a certain threshold
    then the size of the box will be increased.
    ''' 
    if par.step_wise == True:  
        if gen % 2 == 0:
            "Even generations"
            even_list = np.append(even_list,solution)
            
            if gen != 0:
                difference = [np.abs(i[0]-i[1]) for i in zip(even_list,uneven_list)]
                
                if np.mean(np.array(difference)) < 0.5:
                    step+=1
                    #print('Item of boxsize list is now item {}'.format(step))

                even_list = np.zeros([1,12])
                    

        if gen % 2 != 0:
            "Uneven generations"
            uneven_list = np.append(uneven_list,solution)

            difference = [np.abs(i[0]-i[1]) for i in zip(even_list,uneven_list)]
            
            if np.mean(np.array(difference)) < 0.5:
                step+=1
                #print('Item of boxsize list is now item {}'.format(step))
            
            uneven_list = np.zeros([1,12])

    os.chdir(f'{output_dir}')
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
        output.append(f"Average RMSD of parent {solution_idx} in generation {gen}: {population_avg:.4}\n")
        output.append(f"Average RMSD of the lowest RMSDs of parent {solution_idx} in generation {gen}: {population_min_avg:.4}\n")
        output.append("-----------------------------------------------------------------------------------------------------------\n")

        with open('parameter_history', 'a') as f:
            f.write('Parameters generation {}         :  '.format(str(gen))+'  '.join(format(solution[x], ">10.5f") for x in range(0,len(solution)))+'| {:>10.5f}  {:>10.5f}     {:>10.5f}\n'.format(fitness,population_avg, population_min_avg))
    
    else:
        avg_output = np.mean(np.array((avg_list)))
        min_avg_rmsd = np.mean(np.array(min_list))

        fitness = 1 / np.abs(avg_output - desired_output)

        output = [f"-------------------------------------------     PARENT  {solution_idx}      ---------------------------------------------\n"]
        output.append(f"Average RMSD of parent {solution_idx} in generation {gen}: {avg_output:.4}\n")
        output.append(f"Average RMSD of the lowest RMSDs of parent {solution_idx} in generation {gen}: {min_avg_rmsd:.4}\n")
        output.append("-----------------------------------------------------------------------------------------------------------\n")

        with open('parameter_history', 'a') as f:
            f.write('Parameters generation {}         :  '.format(str(gen))+'  '.join(format(solution[x], ">10.5f") for x in range(0,len(solution)))+'| {:>10.5f}  {:>10.5f}     {:>10.5f}\n'.format(fitness, avg_output, min_avg_rmsd))

    gen+=1
    shutil.rmtree(f'{tmp_dir_GA}',ignore_errors=True)

    print(''.join(output))

    return fitness
    

def fitness_wrapper(solution, solution_idx):
    return fitness_func(solution, solution_idx)


class PooledGA(pygad.GA):

    def cal_pop_fitness(self):
        global pool

        pop_fitness = pool.starmap(fitness_wrapper, zip(self.population, list(range(0,len(self.population)))))
        pop_fitness = np.array(pop_fitness)
        return pop_fitness

def train_GA(input_file):
    global pool
    
    global par
    global step
    global gen
    global box_size

    global dir_list
    global tmp_dir

    global even_list
    global uneven_list

    par = Parser(input_file)

    # Global variables
    step = 0
    gen = 0

    ###### Generate Output Dir #######
    if os.path.isdir('output') == False:
        os.mkdir('output')
        os.chdir('output')
    else:
        os.chdir('output')

    if os.path.exists('parameter_history'):
        os.remove('parameter_history')

    if os.path.isdir(f'{output_dir}/first_gen'):
        shutil.rmtree(f'{output_dir}/first_gen', ignore_errors=True)
        os.mkdir(f'{output_dir}/first_gen')
    else:
        os.mkdir(f'{output_dir}/first_gen')

    if os.path.isdir(f'{output_dir}/last_gen'):
        shutil.rmtree(f'{output_dir}/last_gen', ignore_errors=True)
        os.mkdir(f'{output_dir}/last_gen')
    else:
        os.mkdir(f'{output_dir}/last_gen')

    with open('parameter_history', 'a') as f:
        f.write("All old solutions are           :     r_OA        e_OA        r_SA        e_SA        r_HD        e_HD        r_NA        e_NA        r_N         e_N         r_"+par.metal_symbol+"_HD     e_"+par.metal_symbol+"_HD|    fitness    RMSD_AVG   RMSD_MIN_AVG\n")

    # if os.path.isdir('data_set') == False:
    #     # subprocess.call([f'mv {input_dir}/data_set .'], shell=True)
    #     os.chdir('data_set')
    # else:
    os.chdir(f'{input_dir}/data_set')

    # Make list of the protein numbers to iterate over
    dir_list = os.listdir(os.getcwd())
    dir_list = [str(i).replace('protein_','') for i in dir_list]
    dir_list = [int(i) for i in dir_list if convertible(i)]
    dir_list = sorted(dir_list)

    fitness_function = fitness_func

    # Scale factor list
    box_size = np.arange(1.3, 3.1, 0.1)

    num_generations = par.ga_num_generations
    num_parents_mating = par.ga_num_mating

    # initial_population = None

    # sol_per_pop = par.sol_per_pop
    # num_genes = 12
    # gene_space=[{'low': 1, 'high': 3},{'low': 0, 'high': 25},
    #             {'low': 1, 'high': 3},{'low': 0, 'high': 25},
    #             {'low': 1, 'high': 3},{'low': 0, 'high': 25},
    #             {'low': 1, 'high': 3},{'low': 0, 'high': 25},
    #             {'low': 1, 'high': 3},{'low': 0, 'high': 25},
    #             {'low': 1, 'high': 3},{'low': 0, 'high': 25}]

    # parent_selection_type = par.par_type
    # keep_parents = par.keep_par
    # k_tournament = par.k_tour

    # crossover_type = par.crossover_type
    # crossover_probability = par.cross_prob

    # mutation_probability = par.mut_prob
    # mutation_type = mutation_func

    # # Create Class
    # ga_instance = PooledGA(num_generations=num_generations,
    #                        num_parents_mating=num_parents_mating,
    #                        fitness_func=fitness_function,
    #                        sol_per_pop=sol_per_pop,
    #                        num_genes=num_genes,
    #                        gene_space=gene_space,
    #                        parent_selection_type=parent_selection_type,
    #                        keep_parents=keep_parents,
    #                        K_tournament=k_tournament,
    #                        crossover_type=crossover_type,
    #                        crossover_probability=crossover_probability,
    #                        mutation_type=mutation_func,
    #                        save_solutions=True)

    if os.path.isdir(f'{output_dir}/tmp'):
        shutil.rmtree(f'{output_dir}/tmp', ignore_errors=True)
        os.mkdir(f'{output_dir}/tmp')
        os.chdir(f'{output_dir}/tmp')
    else:
        os.mkdir(f'{output_dir}/tmp')
        os.chdir(f'{output_dir}/tmp')

    tmp_dir=os.getcwd()

    # even_list = np.zeros([1,12])
    # uneven_list = np.zeros([1,12])

    # with Pool(processes=sol_per_pop) as pool:
    #     ga_instance.run()

    #     solution, solution_fitness, solution_idx = ga_instance.best_solution()
    #     print("Parameters of the best solution : {solution}\n".format(solution=solution))
    #     print("Fitness value of the best solution = {solution_fitness}\n".format(solution_fitness=solution_fitness))

    # shutil.rmtree(f'{output_dir}/tmp',ignore_errors=True)

    # if os.path.isdir(f'{output_dir}/figures'):
    #     shutil.rmtree(f'{output_dir}/figures', ignore_errors=True)
    #     os.mkdir(f'{output_dir}/figures')
    # else:
    #     os.mkdir(f'{output_dir}/figures')

    # os.chdir(f'{output_dir}/figures')

    # fig.plot_each_protein(dir_list, output_dir)
    # fig.plot_total_conformations(output_dir)
    # fig.plot_parameters(par.metal_symbol, output_dir)

    # print("TRAINING GA COMPLETED")

    parameter_set = [ 2.30177, 0.21898, 1.28486, 16.28000, 1.00407, 5.22678, 1.83437, 23.93663, 1.55722, 8.68594, 1.49070, 11.85050]

    parameter_set = opt.gradient_descent(tmp_dir, dir_list, par, parameter_set, learning_rate=0.5, n_iter=50)

    print('FINAL PARAMETER SET IS: {}'.format(parameter_set))
    return 
