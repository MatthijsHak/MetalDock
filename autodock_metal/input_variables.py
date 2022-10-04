import os
import argparse
import sys, glob
import subprocess
from argparse import RawTextHelpFormatter
import variable_class as vc



def insert_arguments():
    global var
    ###### Argparse ######
    parser = argparse.ArgumentParser(
        prog='MetPar',
        formatter_class=RawTextHelpFormatter,
        description='')


    parser.add_argument("-p", "--parameter_file", type=str, help="Parameter file", metavar="", required=True)
    parser.add_argument("-std", "--standard", action='store_true', help='Dock with standard parameters')
    parser.add_argument("-rmsd", "--rmsd", action='store_true', help='Dock with standard parameters')

    ###### Protein File ######
    parser.add_argument("-pdb","--pdb_file_protein", type=str, help="The pdb file of the protein",required=True)
    parser.add_argument("-pH","--pH",type=str, help='pH of system',required=True)

    ###### Ligand Settings ######
    parser.add_argument("-xyz","--xyz_file_ligand", type=str, help="The xyz file of the ligand",required=True)
    parser.add_argument("-m", "--metal_symbol", type=str, help="Metal atom", required=True)
    parser.add_argument("-c", "--charge_ligand", type=str, help="Charge_ligand", required=True)
    parser.add_argument("-s", "--spin_ligand", type=str, help="Number of unpaired electrons")

    ###### Input Single Point Calculation ######
    parser.add_argument("-basis", "--basis_set", type=str, help="The basis set used in single point", required=True, metavar="")

    functional = parser.add_mutually_exclusive_group(required=True)
    functional.add_argument("-gga", "--gga_functional", type=str, help="Specify the gga functional used", metavar="")
    functional.add_argument("-hybrid", "--hybrid_functional", type=str, help="Specify the hybrid functional", metavar="")

    parser.add_argument("-dispersion", "--dispersion_correction", type=str, help="Specify the dispersion correction", metavar="")

    ###### Docking  Settings######
    parser.add_argument("-random_pos", "--random_position", type=bool, help="The probability of selecting a gene for applying the mutation operation")

    # Box size # 
    box = parser.add_mutually_exclusive_group(required=True)
    box.add_argument("-npts", "--box_size", type=str, help="Manually insert box size")
    box.add_argument("-scale", "--scale_factor", type=float, help="Scale to size of compound.")

    # Docking site # 
    parser.add_argument("-dock_x", "--dock_x", type=str, help="The x-coordinate for docking")
    parser.add_argument("-dock_y", "--dock_y", type=str, help="The y-coordinate for docking")
    parser.add_argument("-dock_z", "--dock_z", type=str, help="The z-coordinate for docking")

    # Docking Algorithm #  
    subparsers = parser.add_subparsers(help='help for subcommand', dest="cmd")
    # GA Settings Docking # 
    parser_GA = subparsers.add_parser("GA_dock", help='Use a genetic algorithm to dock the ligands.')
    parser_GA.add_argument("-ga_pop_size","--ga_pop_size", type=str, help='The population size of the genetic algorithm for the docking procedure.')
    parser_GA.add_argument("-ga_num_evals","--ga_num_evals", type=str, help='The maximum number of energy evluations of the genetic algorithm for the docking procedure.')
    parser_GA.add_argument("-ga_num_generations","--ga_num_generations", type=str, help='The maximum number of generations of the genetic algorithm for the docking procedure.')
    parser_GA.add_argument("-ga_elitism","--ga_elitism", type=str, help='The number of parents kept for the next generation of the genetic algorithm for the docking procedure.')
    parser_GA.add_argument("-ga_mutation_rate","--ga_mutation_rate", type=str, help='The mutation rate of the genetic algorithm for the docking procedure.')
    parser_GA.add_argument("-ga_crossover_rate","--ga_crossover_rate", type=str, help='The crossover rate of the genetic algorithm for the docking procedure.')
    parser_GA.add_argument("-ga_window_size","--ga_window_size", type=str, help='The window size of the genetic algorithm for the docking procedure.')

    # SA Settings Docking # 
    parser_SA = subparsers.add_parser("SA_dock", help='Use the simulated annealing algorithm to dock the ligands.')
    parser_SA.add_argument("-tstep", "--time_step", type=str, help="Type of the mutation operation")
    parser_SA.add_argument("-emax", "--max_initial_energy", type=str, help="The rate at which the mutation changes the gene.")
    parser_SA.add_argument("-rt", "--initial_annealing_temperature", type=str, help="The probability of selecting a gene for applying the mutation operation")
    parser_SA.add_argument("-rtrf", "--temp_reduction_factor", type=str, help="Percentage of genes to mutate")
    parser_SA.add_argument("-runs", "--number_of_runs", type=str, help="Type of the mutation operation")
    parser_SA.add_argument("-cycles", "--max_cycles", type=str, help="The rate at which the mutation changes the gene.")

    ###### Parameter Settings ######
    parser.add_argument("-r_OA", "--r_OA", type=float, help="r_OA")
    parser.add_argument("-e_OA", "--e_OA", type=float, help="e_OA")

    parser.add_argument("-r_SA", "--r_SA", type=float, help="r_SA")
    parser.add_argument("-e_SA", "--e_SA", type=float, help="e_SA")

    parser.add_argument("-r_HD", "--r_HD", type=float, help="r_HD")
    parser.add_argument("-e_HD", "--e_HD", type=float, help="e_HD")

    parser.add_argument("-r_NA", "--r_NA", type=float, help="r_NA")
    parser.add_argument("-e_NA", "--e_NA", type=float, help="e_NA")

    parser.add_argument("-r_N", "--r_N", type=float, help="r_N")
    parser.add_argument("-e_N", "--e_N", type=float, help="e_N")

    parser.add_argument("-r_M", "--r_M", type=float, help="r_M")
    parser.add_argument("-e_M", "--e_M", type=float, help="e_M")

    args = parser.parse_args()

    name_protein = args.pdb_file_protein
    name_protein = name_protein[:-4]

    name_ligand = args.xyz_file_ligand
    name_ligand = name_ligand[:-4]

    metal_cap = args.metal_symbol.upper()

    if args.r_OA and args.e_OA and args.r_SA and args.e_SA and args.r_HD and args.e_HD and args.r_NA and args.e_NA and args.r_N and args.e_N and args.r_M and args.e_M != None:
        args.standard == False 
        print('NON STANDARD PARAMETERS USED!')

    # Random Position #
    if args.random_position == None:
        raise ValueError('You must specifiy if the docking procedure will be performed with an initial random position.')


    if args.cmd == "SA_dock":
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

        docking_simulated_annealing = True
        docking_genetic_algorithm = False
        var = vc.lig_par_dock_SA(args.standard, args.dock_x, args.dock_y, args.dock_z, args.pdb_file_protein, name_protein, args.pH, args.xyz_file_ligand, name_ligand, args.parameter_file, args.scale_factor, args.box_size, args.metal_symbol, 
                                metal_cap, args.charge_ligand, args.spin_ligand, args.basis_set, args.gga_functional, args.hybrid_functional, args.dispersion_correction, docking_simulated_annealing, docking_genetic_algorithm, args.random_position, args.time_step, args.max_initial_energy, 
                                    args.initial_annealing_temperature,args.temp_reduction_factor, args.number_of_runs, args.max_cycles, args.r_OA, args.e_OA, args.r_SA, args.e_SA, args.r_HD, args.e_HD,args.r_NA, args.e_NA,args.r_N, args.e_N, args.r_M, args.e_M)

    if args.cmd == "GA_dock":
        docking_simulated_annealing = False
        docking_genetic_algorithm = True

        var = vc.lig_par_dock_GA(args.standard, args.rmsd, args.dock_x, args.dock_y, args.dock_z, args.pdb_file_protein, name_protein, args.pH, args.xyz_file_ligand, name_ligand, args.parameter_file, args.scale_factor, args.box_size, args.metal_symbol, 
                                    metal_cap, args.charge_ligand, args.spin_ligand, args.basis_set, args.gga_functional, args.hybrid_functional, args.dispersion_correction, docking_simulated_annealing, docking_genetic_algorithm, args.random_position, args.ga_pop_size, args.ga_num_evals, args.ga_num_generations, 
                                        args.ga_elitism, args.ga_mutation_rate, args.ga_crossover_rate, args.ga_window_size, args.r_OA, args.e_OA, args.r_SA, args.e_SA, args.r_HD, args.e_HD,args.r_NA, args.e_NA, args.r_N, args.e_N, args.r_M, args.e_M)
