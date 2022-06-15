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

    parser.add_argument("-r_NA", "--r_nitrogen_A", type=str, help="r_nitrogen_A", metavar="", required=True)
    parser.add_argument("-eps_NA", "--eps_nitrogen_A", type=str, help="eps_nitrogen_A", metavar="", required=True)
    parser.add_argument("-r_N", "--r_nitrogen", type=str, help="r_nitrogen", metavar="", required=True)
    parser.add_argument("-eps_N", "--eps_nitrogen", type=str, help="eps_nitrogen", metavar="", required=True)

    args = parser.parse_args()

    name_protein = args.pdb_file_protein
    name_protein = name_protein.removesuffix('.pdb')

    name_ligand = args.xyz_file_ligand
    name_ligand = name_ligand.removesuffix('.xyz')

    metal_cap = args.metal_symbol.upper()

    global var
    var = vc.lig_par_dock(args.reference_docking, args.dock_x, args.dock_y, args.dock_z, args.rmsd, args.pdb_file_protein, name_protein, args.xyz_file_ligand, name_ligand, args.parameter_file, args.box_size, args.metal_symbol, metal_cap, args.charge_ligand, args.spin_ligand, args.basis_set, args.gga_functional, args.hybrid_functional, args.dispersion_correction, args.r_nitrogen_A, args.eps_nitrogen_A, args.r_nitrogen, args.eps_nitrogen)


