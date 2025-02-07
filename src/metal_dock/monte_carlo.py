import glob
import math
import os
import random
import shutil
import subprocess
import sys
from distutils.dir_util import copy_tree
from itertools import repeat
from multiprocessing import Pool

import numpy as np
import scipy as sc

from src.metal_dock.docking import DockingMC
from src.metal_dock.logger import MetalDockLogger
from src.metal_dock.parser_metal_dock import MCParser
from src.metal_dock.environment_variables import * 

class MonteCarloOptimization:
    def __init__(self, par, input_file):
        self.par = par
        self.input_file = input_file
        self.tmp_dir = None
        self.dir_list = self._prepare_dir_list()
        self.logger = MetalDockLogger()

    def _prepare_dir_list(self):
        """
        Prepare the directory list for the Monte Carlo optimization.
        """
        data_set_dir = self.par.input_dir / 'data_set'
        dir_list = os.listdir(data_set_dir)
        dir_list = [str(i).replace('comopund_', '') for i in dir_list]
        return sorted(dir_list)

    @staticmethod
    def random_sample_continuous():
        """
        Randomly sample a continuous value between 0 and 7.
        """
        return np.random.uniform(low=0, high=7)

    @staticmethod
    def dock_pool(index: int, 
                  input_file: str, 
                  input_dir: str, 
                  tmp_path: str, 
                  dir_list: list, 
                  parameter_set: list):
        """
        Dock the compound in the directory list.

        Args:
            index (int): The index of the compound in the directory list.
            input_file (str): The input file for the Monte Carlo optimization.
            input_dir (str): The input directory for the Monte Carlo optimization.
            tmp_path (str): The temporary path for the Monte Carlo optimization.
            dir_list (list): The directory list for the Monte Carlo optimization.
            parameter_set (list): The parameter set for the Monte Carlo optimization.

        Returns:
            float: The average RMSD of the docked compound.
            list: The print list of the docked compound.
        """
        compound_name = dir_list[index]
        protein_name = f'protein_{index+1}'
        dock_dir_in_path = input_dir / 'data_set' / f'{compound_name}'
        dock_dir_out_path = tmp_path / f'{compound_name}'

        copy_tree(dock_dir_in_path, dock_dir_out_path)

        par = MCParser(config_file=input_file)
        par.method = 'mc'
        par.name_ligand = compound_name
        par.name_protein = protein_name
        par.parameter_set = parameter_set
        par.output_dir = tmp_path
        par.internal_param = False


        docking = DockingMC(par, 
                        metal_complex=None,
                        protein=None,
                        xyz_path=dock_dir_out_path / f'{compound_name}_c.xyz')

        print_list, rmsd_list, avg_rmsd, stdv_rmsd, var_rmsd = docking.run_mc(dock_dir_out_path)

        ##### Fitness function ######
        print_list.insert(0, f'#------------------------------------#')
        print_list.insert(1, f'{compound_name}')
        print_list.insert(2, f'#------------------------------------#')

        return avg_rmsd, print_list

    def _perform_docking_with_parameters(self, 
                                          parameter_set: list):
        """
        Perform the docking with the given parameter set.

        Args:
            parameter_set (list): The parameter set for the Monte Carlo optimization.

        Returns:
            float: The average RMSD of the docked compound.
        """
        with Pool(processes=len(self.dir_list)) as pool:
            results = pool.starmap(MonteCarloOptimization.dock_pool, zip(
                range(len(self.dir_list)), 
                repeat(self.input_file),
                repeat(self.par.input_dir),
                repeat(self.tmp_path),
                repeat(self.dir_list),
                repeat(parameter_set)
            ))

        rmsd_avg_list, print_lists = zip(*results)

        for print_list in print_lists:
            for line in print_list:
                self.logger.info(line)

        rmsd_avg = np.mean(np.array((rmsd_avg_list)))
        return rmsd_avg

    def optimize(self):
        """
        Optimize the parameters with the Monte Carlo optimization scheme.
        """
        self.logger.info('#==============================================================================#')
        self.logger.info("STARTING MONTE CARLO OPTIMIZATION")
        self.logger.info('#==============================================================================#')

        ###### Calculate Box Size #######
        parameter_history = self.par.output_dir / 'parameter_history'
        if parameter_history.exists():
            parameter_history.unlink()

        with open(parameter_history, 'a') as f:
            f.write(f"PARAMETERS         :         e_NA        e_OA        e_SA        e_HD |       RMSD \n")

        self.tmp_path = self.par.output_dir / 'tmp'
        if self.tmp_path.exists():
            shutil.rmtree(self.tmp_path, ignore_errors=True)
        self.tmp_path.mkdir(exist_ok=True)

        new_parameter_set = [2, 2, 2, 2]
        best_parameter_set = [2, 2, 2, 2]
        best_rmsd = self._perform_docking_with_parameters(best_parameter_set)

        i = 0
        step = 1
        while i < self.par.mc_steps:
            for idx in range(len(best_parameter_set)):
                self.logger.info(f"#==============================================================================#")
                self.logger.info(f" MONTE CARLO STEP {step:>4}")
                self.logger.info(f"#==============================================================================#")

                random_eps = self.random_sample_continuous()
                new_parameter_set[idx] = random_eps

                rmsd = self._perform_docking_with_parameters(new_parameter_set)

                self.logger.info(f"RMSD: {rmsd:.5f}")
                self.logger.info(f"PARAMETER SET: {new_parameter_set}\n")

                self.logger.info(f"BEST RMSD: {best_rmsd:.5f}")
                self.logger.info(f"BEST PARAMETER SET: {best_parameter_set}\n")

                if rmsd < best_rmsd:
                    best_parameter_set[idx] = new_parameter_set[idx]
                    best_rmsd = rmsd
                    self.logger.info("PARAMETER SET ACCEPTED")
                    self.logger.info(f"NEW BEST PARAMETER SET: {best_parameter_set}\n")
                    self._log_parameter_history(step, best_parameter_set, best_rmsd)
                else:
                    diff_RMSD = best_rmsd - rmsd
                    acceptance_probability = np.exp(diff_RMSD)
                    random_number = random.uniform(0, 1)

                    if random_number < acceptance_probability:
                        best_parameter_set[idx] = new_parameter_set[idx]
                        best_rmsd = rmsd
                        self.logger.info("ACCEPTED WITH HIGHER RMSD")
                        self.logger.info(f"NEW PARAMETER SET: {best_parameter_set}\n")
                        self._log_parameter_history(step, best_parameter_set, best_rmsd)
                    else:
                        new_parameter_set[idx] = best_parameter_set[idx]
                        self.logger.info("PARAMETER SET DENIED")

                step += 1
            i += 1

        self._finalize_optimization()

    def _log_parameter_history(self, 
                              step: int, 
                              best_parameter_set: list, 
                              best_rmsd: float):
        """
        Log the parameter history.

        Args:
            step (int): The step of the Monte Carlo optimization.
            best_parameter_set (list): The best parameter set.
            best_rmsd (float): The best RMSD.
        """
        with open(self.par.output_dir / 'parameter_history', 'a') as f:
            f.write(f"STEP {step:>6}        :   " + "  ".join(f"{param:>10.5f}" for param in best_parameter_set) + f" | {best_rmsd:>10.5f}\n")

    def _finalize_optimization(self):
        """
        Finalize the optimization.
        """
        with open(self.par.output_dir / 'parameter_history', 'r') as fin:
            lines = [line.strip().split() for line in fin]
            rmsd = [float(line[-1]) for line in lines[1:]]
            min_idx = np.argmin(rmsd) + 1
            line = lines[min_idx]

        shutil.rmtree(f'{self.tmp_path}')
        os.chdir(f'{self.par.output_dir}')
        self.logger.info(f"#==============================================================================#")
        self.logger.info(f"BEST RMSD: {line[-1]}")
        self.logger.info(f"BEST PARAMETERS: e_NA {line[3]} kcal/mol; e_OA {line[4]} kcal/mol; e_SA {line[5]} kcal/mol; e_HD {line[6]} kcal/mol")
        self.logger.info(f"#==============================================================================#")