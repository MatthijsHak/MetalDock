import os
import sys
import shutil
import itertools
import numpy as np

from ase.io import read
from ase.calculators.gaussian import Gaussian, GaussianOptimizer
from src.metal_dock.logger import MetalDockLogger

class GaussianEngine:
    def __init__(self, par, metal_complex):
        self.par = par
        self.metal_complex = metal_complex
        self.logger = MetalDockLogger()

    def run(self):
        """
        Run the Gaussian engine.
        """
        self.logger.info('\nGAUSSIAN ENGINE CURRENTLY HAS A BUG DUE TO ASE MODULE')
        self.logger.info("THERE WILL BE AN ERRROR AFTER RUNNING THE QM CALCULATION")
        self.logger.info("PLEASE VERIFY IF THE CALCULATION WAS SUCCESSFUL BY CHECKING THE LOG FILE\n")
        self.logger.info("IF THE CALCULATION WAS SUCCESSFUL, PLEASE RUN METALDOCK AGAIN")
        self.logger.info("EVERYTHING SHOULD THEN STILL RUN SMOOTHLY\n")

        if self.par.geom_opt:
            self._run_geom_opt()
        else:
            self._run_single_point()

        return os.getcwd(), self.energy

    def _run_geom_opt(self):
        """
        Run the geometry optimization.
        """
        os.makedirs(self.par.output_dir / 'QM' / 'geom_opt', exist_ok=True)

        shutil.copy(self.par.xyz_file, self.par.output_dir / 'QM' / 'geom_opt' / self.par.xyz_file)

        log_file = self.par.output_dir / 'QM' / 'geom_opt' / 'geom_opt.log'
        chk_file = self.par.output_dir / 'QM' / 'geom_opt' / 'geom_opt.chk'
        if not os.path.exists(chk_file):
            self._gaussian_geom_opt()

        if self._gaussian_opt_converged(log_file):
            xyz_file = self.par.output_dir / 'QM' / 'geom_opt' / 'output.xyz'
            cm5_path = self.par.output_dir / 'QM' / 'geom_opt' / 'CM5_charges'
            self.energy = self._gaussian_extract_energy(log_file)
            self.metal_complex.energy = self.energy

            self._gaussian_extract_CM5(log_file, xyz_file, cm5_path)
            self._gaussian_extract_charges('geom_opt')

        os.chdir(self.par.output_dir)
        self.metal_complex.create_mol_graph('geom_opt')
        self.metal_complex.add_charges_to_graph()
        self._gaussian_extract_bond_orders(log_file, xyz_file)

    def _run_single_point(self):
        """
        Run the single point calculation.
        """
        os.makedirs(self.par.output_dir / 'QM' / 'single_point', exist_ok=True)

        shutil.copy(self.par.xyz_file, self.par.output_dir / 'QM' / 'single_point' / self.par.xyz_file)

        chk_file = self.par.output_dir / 'QM' / 'single_point' / 'single_point.chk'
        log_file = self.par.output_dir / 'QM' / 'single_point' / 'single_point.log'
        if not os.path.exists(chk_file):
            self._gaussian_sp()

        if self._gaussian_sp_converged(log_file):
            xyz_file = self.par.output_dir / 'QM' / 'single_point' / 'output.xyz'
            cm5_path = self.par.output_dir / 'QM' / 'single_point' / 'CM5_charges'
            
            self.energy = self._gaussian_extract_energy(log_file)
            self.metal_complex.energy = self.energy

            self._gaussian_extract_CM5(log_file, xyz_file, cm5_path)
            self._gaussian_extract_charges('single_point')

        os.chdir(self.par.output_dir)
        self.metal_complex.create_mol_graph('single_point')
        self.metal_complex.add_charges_to_graph()
        self._gaussian_extract_bond_orders(log_file, xyz_file)

    def _gaussian_extract_energy(self, log_file):
        """
        Extract the energy from the log file.
        """
        with open(log_file, 'r') as fin:
            for line in fin:
                if line.startswith(' SCF Done:'):
                    energy = line.split()[4]
                    return energy

    def _gaussian_extract_CM5(self, log_file: str, xyz_file: str, cm5_path: str):
        """
        Extract the CM5 charges from the log file.

        Args:
            log_file (str): The path to the log file.
            xyz_file (str): The path to the xyz file.
            cm5_path (str): The path to the cm5 file.
        """
        mol = read(xyz_file)
        N = len(mol.positions)

        with open(log_file) as log:
            for line in log:
                if ' Hirshfeld charges, spin densities, dipoles, and CM5 charges' in line:
                    fin_lines = list(itertools.islice(log, N + 1))
                    fin_lines = [line.strip().split() for line in fin_lines]
                    with open(cm5_path, 'w') as fout:
                        for i in fin_lines[1:]:
                            fout.write('{} {}\n'.format(i[1], i[7]))

    def _gaussian_extract_charges(self, run_type):
        """
        Extract the CM5 charges from the log file.
        """
        charges = {}
        atom_id = 0
        with open(self.par.output_dir / 'QM' / run_type / 'CM5_charges', 'r') as file:
            for line in file:
                line = line.strip()
                if line:  # Ensure the line is not empty
                    atom_id += 1
                    parts = line.split()
                    atom = parts[0].upper()  # Atom type
                    charge = float(parts[1])  # Charge value
                    charges[f'{atom}{atom_id}'] = charge

        self.metal_complex.charges = charges

    def _gaussian_extract_bond_orders(self, log_file: str, xyz_file: str):
        """
        Extract the bond orders from the log file.

        Args:
            log_file (str): The path to the log file.
            xyz_file (str): The path to the xyz file.
        """
        mol = read(xyz_file)
        n_atoms = len(mol.positions)

        with open(log_file, 'r') as file:
            lines = file.readlines()
        
        # Initialize variables
        extract_lines = []
        recording = False
        
        for line in lines:
            if ' Atomic Valencies and Mayer Atomic Bond Orders:' in line:
                recording = True
                continue
            elif ' Lowdin Atomic Charges:' in line:
                recording = False
                break
            
            if recording:
                extract_lines.append(line)

        # Remove any empty lines and strip whitespace
        lines = [line.strip() for line in extract_lines if line.strip()]

        matrix = np.zeros((n_atoms, n_atoms))
        n_chunks = n_atoms // 6
        rest_columns = n_atoms % 6

        for chunk in range(n_chunks + 1):  # +1 to include the final smaller chunk
            start_line = 1 + chunk * (n_atoms+1)  # Skip header line and previous chunks
            end_line = start_line + n_atoms  # Each chunk has 74 data lines
            
            if chunk < n_chunks:
                columns = 6
            else:
                columns = rest_columns

            for col in range(columns):
                for row in range(n_atoms):
                    line = lines[start_line + row].split()
                    value = float(line[col + 2])  # +2 to skip index and atom symbol
                    matrix[row][chunk * 6 + col] = value

        # Create bond orders that are larger than 0.2 in the following format (atom1, atom2): bond_order
        bond_orders = {}
        for i in range(n_atoms):
            for j in range(i+1, n_atoms):
                if matrix[i][j] > 0.2:  # Only consider bond orders greater than 0.2
                    bond_orders[(i, j)] = matrix[i][j]

        # Update the metal_complex graph with bond orders
        for (atom1, atom2), order in bond_orders.items():
            if self.metal_complex.graph.has_edge(atom1, atom2):
                self.metal_complex.graph[atom1][atom2]['bond_order'] = order

        # Remove edges without bond orders from the graph
        for edge in list(self.metal_complex.graph.edges()):
            sorted_edge = tuple(sorted(edge))
            if sorted_edge not in bond_orders:
                self.metal_complex.graph.remove_edge(*edge)

    def _gaussian_opt_converged(self, log_file: str):
        """
        Check if the geometry optimization has converged.

        Args:
            log_file (str): The path to the log file.

        Returns:
            bool: True if the geometry optimization has converged, False otherwise.
        """
        with open(log_file) as log:
            if 'Optimization completed.' in log.read():
                self.logger.info('GEOMETRY CONVERGED')
                return True
            else:
                self.logger.info('GEOMETRY NOT CONVERGED\n CHECK .log FILE IN QM DIRECTORY\nDELETE .chk, .com & .log FILES TO RERUN\n')
                sys.exit()

    def _gaussian_sp_converged(self, log_file: str):
        """
        Check if the single point calculation has converged.

        Args:
            log_file (str): The path to the log file.

        Returns:
            bool: True if the single point calculation has converged, False otherwise.
        """
        with open(log_file) as log:
            if 'SCF Done' in log.read():
                self.logger.info('\nSINGLE POINT SUCCESSFULLY PERFORMED\n')
                return True
            else:
                self.logger.info('\nSINGLE POINT NOT SUCCESSFUL\nCHECK .log FILE IN QM DIRECTORY\nDELETE .chk, .com & .log FILES TO RERUN\n')
                sys.exit()

    def _gaussian_geom_opt(self):
        """
        Run the geometry optimization.
        """
        os.chdir(self.par.output_dir / 'QM' / 'geom_opt')
        M = 2 * (self.par.spin * 0.5) + 1

        mol = read(self.par.xyz_file)

        if self.par.solvent and self.par.dispersion:
            s = Gaussian(label='geom_opt',
                         nprocshared=self.par.ncpu,
                         mem=f'{self.par.memory}MB',
                         chk='geom_opt.chk',
                         xc=self.par.functional,
                         charge=self.par.charge,
                         mult=M,
                         basis=self.par.basis_set,
                         pop='Hirshfeld',
                         SCRF=f'PCM, solvent={self.par.solvent}',
                         EmpiricalDispersion=self.par.dispersion,
                         ioplist=['6/80=1'])

        elif not self.par.solvent and self.par.dispersion:
            s = Gaussian(label='geom_opt',
                         nprocshared=self.par.ncpu,
                         mem=f'{self.par.memory}MB',
                         chk='geom_opt.chk',
                         xc=self.par.functional,
                         charge=self.par.charge,
                         mult=M,
                         basis=self.par.basis_set,
                         pop='Hirshfeld',
                         EmpiricalDispersion=self.par.dispersion,
                         ioplist=['6/80=1'])

        elif self.par.solvent and not self.par.dispersion:
            s = Gaussian(label='geom_opt',
                         nprocshared=self.par.ncpu,
                         mem=f'{self.par.memory}MB',
                         chk='geom_opt.chk',
                         xc=self.par.functional,
                         charge=self.par.charge,
                         mult=M,
                         basis=self.par.basis_set,
                         pop='Hirshfeld',
                         SCRF=f'PCM, solvent={self.par.solvent}')

        else:
            s = Gaussian(label='geom_opt',
                         nprocshared=self.par.ncpu,
                         mem=f'{self.par.memory}MB',
                         chk='geom_opt.chk',
                         xc=self.par.functional,
                         charge=self.par.charge,
                         mult=M,
                         basis=self.par.basis_set,
                         pop='Hirshfeld')

        opt = GaussianOptimizer(mol, s)
        opt.run(fmax='tight')
        mol.write(self.par.output_dir / 'QM' / 'geom_opt' / 'output.xyz')

    def _gaussian_sp(self):
        """
        Run the single point calculation.
        """
        os.chdir(self.par.output_dir / 'QM' / 'single_point')
        M = 2 * (self.par.spin * 0.5) + 1

        mol = read(self.par.xyz_file)
        if self.par.solvent and self.par.dispersion:
            s = Gaussian(label='single_point',
                         nprocshared=self.par.ncpu,
                         mem=f'{self.par.memory}MB',
                         chk='single_point.chk',
                         xc=self.par.functional,
                         charge=self.par.charge,
                         mult=M,
                         basis=self.par.basis_set,
                         pop='Hirshfeld',
                         SCRF=f'PCM, solvent={self.par.solvent}',
                         EmpiricalDispersion=self.par.dispersion,
                         ioplist=['6/80=1'])

        elif not self.par.solvent and self.par.dispersion:
            s = Gaussian(label='single_point',
                         nprocshared=self.par.ncpu,
                         mem=f'{self.par.memory}MB',
                         chk='single_point.chk',
                         xc=self.par.functional,
                         charge=self.par.charge,
                         mult=M,
                         basis=self.par.basis_set,
                         pop='Hirshfeld',
                         EmpiricalDispersion=self.par.dispersion,
                         ioplist=['6/80=1'])

        elif self.par.solvent and not self.par.dispersion:
            s = Gaussian(label='single_point',
                         nprocshared=self.par.ncpu,
                         mem=f'{self.par.memory}MB',
                         chk='single_point.chk',
                         xc=self.par.functional,
                         charge=self.par.charge,
                         mult=M,
                         basis=self.par.basis_set,
                         pop='Hirshfeld',
                         SCRF=f'PCM, solvent={self.par.solvent}',
                         ioplist=['6/80=1'])

        else:
            s = Gaussian(label='single_point',
                         nprocshared=self.par.ncpu,
                         mem=f'{self.par.memory}MB',
                         chk='single_point.chk',
                         xc=self.par.functional,
                         charge=self.par.charge,
                         mult=M,
                         basis=self.par.basis_set,
                         pop='Hirshfeld',
                         ioplist=['6/80=1'])

        mol.calc = s
        mol.get_potential_energy()
        mol.write(self.par.output_dir / 'QM' / 'single_point' / 'output.xyz')