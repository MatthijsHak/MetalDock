import os
import sys
import argparse

from pathlib import Path

from src.metal_dock.parser_metal_dock import DockParser, MCParser
from src.metal_dock.protein import Protein
from src.metal_dock.metal_complex import MetalComplex
from src.metal_dock.docking import Docking
from src.metal_dock.monte_carlo import MonteCarloOptimization
from src.metal_dock.logger import MetalDockLogger

def docking(par: DockParser):
    """
    Docking function.

    Args:
        par (DockParser): The parser object.
    """
    input_dir = Path.cwd()
    par.input_dir = input_dir
    par.output_dir = input_dir / 'output'
    #=== OUTPUT DIRECTORY ===#
    par.output_dir.mkdir(exist_ok=True)

    #=== FILE PREPARATION ===#
    file_prep_dir = par.output_dir / 'file_prep'
    file_prep_dir.mkdir(exist_ok=True)

    #=== PROTEIN INITIALIZATION ===#
    # initialize the protein
    protein = Protein(par)
    protein.create_pdbqt_file()

    #=== METAL COMPLEX INITIALIZATION ===#
    metal_complex = MetalComplex(par)
    metal_complex.canonicalize_ligand()

    # run the qm engine
    metal_complex.qm_engine.run()

    # create the pdbqt file
    metal_complex.create_ligand_pdbqt_file()

    #=== DOCKING ===#
    docking_dir = par.output_dir / 'docking'
    docking_dir.mkdir(exist_ok=True)

    docking = Docking(par, metal_complex, protein)
    docking.run()

    results_dir = par.output_dir / 'results'
    results_dir.mkdir(exist_ok=True)

    docking.analyze_results()

def optimize_MC(par: MCParser, input_file: str):
    """
    Function to optimize the parameters with the Monte Carlo optimization.

    Args:
        par (MCParser): The parser object.
        input_file (str): The input file.
    """
    input_dir = Path.cwd()
    par.input_dir = input_dir
    par.output_dir = input_dir / 'output'
    #=== OUTPUT DIRECTORY ===#
    par.output_dir.mkdir(exist_ok=True)

    optimize_MC = MonteCarloOptimization(par, input_file)
    optimize_MC.optimize()

def main():
    parser = argparse.ArgumentParser(description="Docking of organometallic compounds")
    parser.add_argument("-i","--input_file", help="Input file of .ini format specifying all parameters", required=True)
    parser.add_argument("-m","--method", choices=['dock', 'mc'], help="Method to use", required=True)
    args = parser.parse_args()

    logger = MetalDockLogger()
    logger.setup_file_logger(Path.cwd())

    logger.info('#==============================================================================#')
    logger.info("STARTING METALDOCK")
    logger.info('#==============================================================================#\n')

    # if input file not found give errror
    if not os.path.exists(args.input_file):
        logger = MetalDockLogger()
        logger.info("Input file not found")
        sys.exit()

    if args.method.lower() == 'dock':
        par = DockParser(args.input_file)
        par.method = 'dock'
        docking(par)
    elif args.method.lower() == 'mc':
        par = MCParser(args.input_file)
        par.method = 'mc'
        optimize_MC(par, args.input_file)
    else:
        logger = MetalDockLogger()
        logger.info("SPECIFY ONE OF THE TWO OPTIONS FOR MetalDock")
        logger.info("(1) dock")
        logger.info("(2) mc")

if __name__== '__main__':
    main()
