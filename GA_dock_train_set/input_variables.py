import os
import argparse
import sys
import subprocess
from argparse import RawTextHelpFormatter
import variable_class as vc



def insert_arguments():
    ###### Argparse ######
    parser = argparse.ArgumentParser(
        prog='MetPar',
        formatter_class=RawTextHelpFormatter,
        description='')


    ###### Ligand Settings ######
    parser.add_argument("-m", "--metal_symbol", type=str, help="Metal atom", required=True)

    ###### Parameter Settings ######
    parser.add_argument("-p", "--parameter_file", type=str, help="Parameter file", required=True)

    # Box size # 
    box = parser.add_mutually_exclusive_group(required=True)
    box.add_argument("-npts", "--box_size", type=str, help="Manually insert box size")
    box.add_argument("-scale", "--scale_factor", type=float, help="Scale to size of compound.")

    docking_method = parser.add_mutually_exclusive_group(required=True)
    docking_method.add_argument('-GA','--docking_genetic_algorithm', action='store_true')
    docking_method.add_argument('-SA','--docking_simulated_annealing', action='store_true')

    ###### GA Settings Metal Parameters######
    parser.add_argument("-num_gen", "--num_generations", type=int, help="Number of generations within population")
    parser.add_argument("-num_mat", "--num_parents_mating", type=int, help="Number of solutions to be selected as parents")
    parser.add_argument("-sol", "--sol_per_pop", type=int, help="Number of solutions within population")

    # Selection #
    parser.add_argument("-par_type", "--parent_selection_type", type=str, help="The parent selection type")
    parser.add_argument("-keep_par", "--keep_parents", type=int, help="Number of parents to keep in the current population")
    parser.add_argument("-k_tour", "--k_tournament", type=int, help="Parameter file")

    # Crossover #
    parser.add_argument("-cross", "--crossover_type", type=str, help="Type of the crossover operation")
    parser.add_argument("-cross_prob", "--crossover_prob", type=float, help="The probability of selecting a parent for applying the crossover operation")

    # Mutation #
    parser.add_argument("-mut_prob", "--mutation_prob", type=float, help="The probability of selecting a gene for applying the mutation operation")

    ###### Docking  Settings######
    parser.add_argument("-random_pos", "--random_position", type=bool, help="The probability of selecting a gene for applying the mutation operation")

    # SA Settings Docking # 
    parser.add_argument("-tstep", "--time_step", type=str, help="Type of the mutation operation")
    parser.add_argument("-emax", "--max_initial_energy", type=str, help="The rate at which the mutation changes the gene.")
    parser.add_argument("-rt", "--initial_annealing_temperature", type=str, help="The probability of selecting a gene for applying the mutation operation")
    parser.add_argument("-rtrf", "--temp_reduction_factor", type=str, help="Percentage of genes to mutate")
    parser.add_argument("-runs", "--number_of_runs", type=str, help="Type of the mutation operation")
    parser.add_argument("-cycles", "--max_cycles", type=str, help="The rate at which the mutation changes the gene.")

    # GA Settings Docking # 

    args = parser.parse_args()

    metal_cap = args.metal_symbol.upper()

    global var


    if args.num_generations == None:
        raise ValueError('You must specifiy the number of generations to use the gentic algorithm as docking method.')

    if args.num_parents_mating == None:
        raise ValueError('You must specifiy the number of parents mating to use the gentic algorithm as docking method.')
        
    if args.sol_per_pop == None:
        raise ValueError('You must specifiy the solutions per population to use the gentic algorithm as docking method.')

    # Selection #
    if args.parent_selection_type == None:
        raise ValueError('You must specifiy the parent selection type to use the gentic algorithm as docking method.')

    if args.parent_selection_type == 'tournament' and args.K_tournament == None:
        raise ValueError('The number of parents that will be used in the tournament selection is not specified.')

    if args.keep_parents == None:
        raise ValueError('You must specifiy the number of parents to keep to use the gentic algorithm as docking method.')

    # Crossover #
    if args.crossover_type == None:
        raise ValueError('You must specifiy a crossover operation to use the gentic algorithm as docking method.')

    if args.crossover_prob == None:
        raise ValueError('You must specifiy the parent selection type to use the gentic algorithm as docking method.')

    # Mutation #
    if args.mutation_prob == None:
        raise ValueError('You must specifiy the mutation probability  to use the gentic algorithm as docking method.')

    # Random Position #
    if args.random_position == None:
        raise ValueError('You must specifiy if the docking procedure will be performed with an initial random position.')

    if args.docking_simulated_annealing == True:
        if args.time_step == None:
            raise ValueError('You must specifiy the time step to use simulated annealing as docking method.')

        if args.max_initial_energy == None:
            raise ValueError('You must specifiy the max initial energy to use simulated annealing as docking method.')

        if args.initial_annealing_temperature == None:
            raise ValueError('You must specifiy the initial annealing temperature to use simulated annealing as docking method.')

        if args.temp_reduction_factor == None:
            raise ValueError('You must specifiy the temparture reduction factor to use simulated annealing as docking method.')

        if args.number_of_runs == None:
            raise ValueError('You must specifiy the number of runs to use simulated annealing as docking method.')

        if args.max_cycles == None:
            raise ValueError('You must specifiy the maximum number of cylces for each run to use simulated annealing as docking method.')
        
    var = vc.variables(args.parameter_file, args.box_size, args.scale_factor, args.metal_symbol, metal_cap, args.num_generations, args.num_parents_mating, args.sol_per_pop, args.parent_selection_type, 
                            args.keep_parents, args.k_tournament, args.crossover_type, args.crossover_prob, args.mutation_prob , args.docking_simulated_annealing, args.docking_genetic_algorithm, 
                                    args.random_position, args.time_step, args.max_initial_energy, args.initial_annealing_temperature, args.temp_reduction_factor, args.number_of_runs, args.max_cycles)
                                         


