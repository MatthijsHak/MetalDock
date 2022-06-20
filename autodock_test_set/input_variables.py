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
    parser.add_argument("-p", "--parameter_file", type=str, help="Parameter file", metavar="", required=True)
    parser.add_argument("-npts", "--box_size", type=str, help="Boxsize", metavar="", required=True)

    parser.add_argument("-r_OA", "--r_OA", type=str, help="r_OA", required=True)
    parser.add_argument("-e_OA", "--e_OA", type=str, help="e_OA", required=True)

    parser.add_argument("-r_SA", "--r_SA", type=str, help="r_OA", required=True)
    parser.add_argument("-e_SA", "--e_SA", type=str, help="e_SA", required=True)

    parser.add_argument("-r_HD", "--r_HD", type=str, help="r_HD", required=True)
    parser.add_argument("-e_HD", "--e_HD", type=str, help="e_HD", required=True)

    parser.add_argument("-r_NA", "--r_NA", type=str, help="r_NA", required=True)
    parser.add_argument("-e_NA", "--e_NA", type=str, help="e_NA", required=True)

    parser.add_argument("-r_N", "--r_N", type=str, help="r_N", required=True)
    parser.add_argument("-e_N", "--e_N", type=str, help="e_N", required=True)

    parser.add_argument("-r_Ru_Ru", "--r_Ru_Ru", type=str, help="r_N", required=True)
    parser.add_argument("-e_Ru_Ru", "--e_Ru_Ru", type=str, help="e_N", required=True)


    args = parser.parse_args()

    metal_cap = args.metal_symbol.upper()

    global var
    var = vc.lig_par_dock(args.parameter_file, args.box_size, args.metal_symbol, metal_cap, args.r_OA, args.e_OA, args.r_SA, args.e_SA,args.r_HD, args.e_HD,args.r_NA, args.e_NA,args.r_N, args.e_N, args.r_Ru_Ru, args.e_Ru_Ru )


