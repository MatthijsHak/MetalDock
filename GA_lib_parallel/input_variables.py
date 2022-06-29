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
    parser.add_argument("-npts", "--box_size", type=str, help="Boxsize", required=True)

    ###### GA Settings ######
    parser.add_argument("-num_gen", "--num_generations", type=int, help="Number of generations within population", required=True)
    parser.add_argument("-num_mat", "--num_parents_mating", type=int, help="Number of solutions to be selected as parents", required=True)

    parser.add_argument("-sol", "--sol_per_pop", type=int, help="Number of solutions within population", required=True)

    # Selection #
    parser.add_argument("-par_type", "--parent_selection_type", type=str, help="The parent selection type", required=True)
    parser.add_argument("-keep_par", "--keep_parents", type=int, help="Number of parents to keep in the current population", required=True)
    parser.add_argument("-K_tour", "--K_tournament", type=int, help="Parameter file")

    # Crossover #
    parser.add_argument("-cross", "--crossover_type", type=str, help="Type of the crossover operation", required=True)
    parser.add_argument("-cross_prob", "--crossover_prob", type=float, help="The probability of selecting a parent for applying the crossover operation")

    # Mutation #
    parser.add_argument("-mut", "--mutation_type", type=float, help="Type of the mutation operation")
    parser.add_argument("-mut_rate", "--mutation_rate", type=float, help="The rate at which the mutation changes the gene.")
    parser.add_argument("-mut_prob", "--mutation_prob", type=float, help="The probability of selecting a gene for applying the mutation operation")
    parser.add_argument("-mut_perc", "--mutation_percent", type=int, help="Percentage of genes to mutate")

    args = parser.parse_args()

    metal_cap = args.metal_symbol.upper()

    global var
    var = vc.lig_par_dock(args.parameter_file, args.box_size, args.metal_symbol, metal_cap, args.num_generations, args.num_parents_mating, args.sol_per_pop, args.parent_selection_type, args.keep_parents, args.K_tournament, args.crossover_type, args.crossover_prob, args.mutation_type, args.mutation_rate, args.mutation_prob, args.mutation_percent)


