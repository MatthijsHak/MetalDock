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

    ###### Docking  Settings######
    parser.add_argument("-random_pos", "--random_position", type=bool, help="The probability of selecting a gene for applying the mutation operation")

    # SA Settings Docking # 
    parser.add_argument("-tstep", "--time_step", type=str, help="Type of the mutation operation")
    parser.add_argument("-emax", "--max_initial_energy", type=str, help="The rate at which the mutation changes the gene.")
    parser.add_argument("-rt", "--initial_annealing_temperature", type=str, help="The probability of selecting a gene for applying the mutation operation")
    parser.add_argument("-rtrf", "--temp_reduction_factor", type=str, help="Percentage of genes to mutate")
    parser.add_argument("-runs", "--number_of_runs", type=str, help="Type of the mutation operation")
    parser.add_argument("-cycles", "--max_cycles", type=str, help="The rate at which the mutation changes the gene.")


    # Optimized parameters from GA #
    parser.add_argument("-r_OA", "--r_OA", type=float, help="r_OA")
    parser.add_argument("-e_OA", "--e_OA", type=float, help="e_OA")

    parser.add_argument("-r_SA", "--r_SA", type=float, help="r_OA")
    parser.add_argument("-e_SA", "--e_SA", type=float, help="e_SA")

    parser.add_argument("-r_HD", "--r_HD", type=float, help="r_HD")
    parser.add_argument("-e_HD", "--e_HD", type=float, help="e_HD")

    parser.add_argument("-r_NA", "--r_NA", type=float, help="r_NA")
    parser.add_argument("-e_NA", "--e_NA", type=float, help="e_NA")

    parser.add_argument("-r_N", "--r_N", type=float, help="r_N")
    parser.add_argument("-e_N", "--e_N", type=float, help="e_N")

    parser.add_argument("-r_M", "--r_M", type=float, help="r_N")
    parser.add_argument("-e_M", "--e_M", type=float, help="e_N")


    args = parser.parse_args()

    metal_cap = args.metal_symbol.upper()


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
        
    global var
    var = vc.lig_par_dock(args.metal_symbol, metal_cap, args.parameter_file, args.box_size, args.scale_factor, args.docking_genetic_algorithm, args.docking_simulated_annealing,
                                args.random_position, args.time_step, args.max_initial_energy, args.initial_annealing_temperature, args.temp_reduction_factor, args.number_of_runs,
                                    args.max_cycles, args.r_OA, args.e_OA, args.r_SA, args.e_SA,args.r_HD, args.e_HD,args.r_NA, args.e_NA,args.r_N, args.e_N, args.r_M, args.e_M)


