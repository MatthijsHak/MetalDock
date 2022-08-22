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


    ##### Docking type #####
    parser.add_argument("-ref", "--reference_docking", type=bool, help="Docking is done to calculate the rmsd of reference structure")

    ##### Docking type #####
    parser.add_argument("-dock_x", "--dock_x", type=str, help="The x-coordinate for docking")
    parser.add_argument("-dock_y", "--dock_y", type=str, help="The y-coordinate for docking")
    parser.add_argument("-dock_z", "--dock_z", type=str, help="The z-coordinate for docking")

    ###### Redock ######
    parser.add_argument("-rmsd", "--rmsd", type=bool, help="Only calculate RMSD of docked ligands")

    ###### Protein File ######
    parser.add_argument("-pdb","--pdb_file_protein", type=str, help="The pdb file of the protein",required=True)
    parser.add_argument("-pH","--pH",type=str, help='pH of system',required=True)

    ###### Ligand Settings ######
    parser.add_argument("-xyz","--xyz_file_ligand", type=str, help="The xyz file of the ligand",required=True)
    parser.add_argument("-m", "--metal_symbol", type=str, help="Metal atom", required=True)
    parser.add_argument("-c", "--charge_ligand", type=str, help="Charge_ligand", required=True)
    parser.add_argument("-s", "--spin_ligand", type=str, help="Number of unpaired electrons")

    ###### Input Single Point & Hessian Calculation ######
    parser.add_argument("-basis", "--basis_set", type=str, help="The basis set used in single point", required=True, metavar="")

    functional = parser.add_mutually_exclusive_group(required=True)
    functional.add_argument("-gga", "--gga_functional", type=str, help="Specify the gga functional used", metavar="")
    functional.add_argument("-hybrid", "--hybrid_functional", type=str, help="Specify the hybrid functional", metavar="")

    parser.add_argument("-dispersion", "--dispersion_correction", type=str, help="Specify the dispersion correction", metavar="")


    ###### Parameter Settings ######
    parser.add_argument("-p", "--parameter_file", type=str, help="Parameter file", metavar="", required=True)
    parser.add_argument("-npts", "--box_size", type=str, help="Boxsize", metavar="", required=True)

    parser.add_argument("-r_OA", "--r_OA", type=str, help="r_OA")
    parser.add_argument("-e_OA", "--e_OA", type=str, help="e_OA")

    parser.add_argument("-r_SA", "--r_SA", type=str, help="r_SA")
    parser.add_argument("-e_SA", "--e_SA", type=str, help="e_SA")

    parser.add_argument("-r_HD", "--r_HD", type=str, help="r_HD")
    parser.add_argument("-e_HD", "--e_HD", type=str, help="e_HD")

    parser.add_argument("-r_NA", "--r_NA", type=str, help="r_NA")
    parser.add_argument("-e_NA", "--e_NA", type=str, help="e_NA")

    parser.add_argument("-r_N", "--r_N", type=str, help="r_N")
    parser.add_argument("-e_N", "--e_N", type=str, help="e_N")

    parser.add_argument("-r_M", "--r_M", type=str, help="r_M")
    parser.add_argument("-e_M", "--e_M", type=str, help="e_M")


    args = parser.parse_args()

    if args.r_OA == None:
        args.r_OA = "2.0"
    if args.e_OA == None:
        args.e_OA = "10.0"
    if args.r_SA == None:
        args.r_SA = "2.0"
    if args.e_SA == None:
        args.e_SA = "10.0"
    if args.r_HD == None:
        args.r_HD = "2.0"
    if args.e_HD == None:
        args.e_HD = "10.0"
    if args.r_NA == None:
        args.r_NA = "2.0"
    if args.e_NA == None:
        args.e_NA = "10.0"
    if args.r_N == None:
        args.r_N = "2.0"
    if args.e_N == None:
        args.e_N = "10.0"
    if args.r_M == None:
        args.r_M = "2.0"
    if args.e_M == None:
        args.e_M = "10.0"

    name_protein = args.pdb_file_protein
    name_protein = name_protein.removesuffix('.pdb')

    name_ligand = args.xyz_file_ligand
    name_ligand = name_ligand.removesuffix('.xyz')

    metal_cap = args.metal_symbol.upper()

    global var
    var = vc.lig_par_dock(args.reference_docking, args.dock_x, args.dock_y, args.dock_z, args.rmsd, args.pdb_file_protein, name_protein, args.pH, args.xyz_file_ligand, name_ligand, args.parameter_file, args.box_size, args.metal_symbol, metal_cap, args.charge_ligand, args.spin_ligand, args.basis_set, args.gga_functional, args.hybrid_functional, args.dispersion_correction, args.r_OA, args.e_OA, args.r_SA, args.e_SA,args.r_HD, args.e_HD,args.r_NA, args.e_NA,args.r_N, args.e_N, args.r_Ru_Ru, args.e_Ru_Ru)

