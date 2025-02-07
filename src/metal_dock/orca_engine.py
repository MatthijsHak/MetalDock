import os
import sys
import shutil
import re
import json

import numpy as np
import pandas as pd

from ase.io import read
from ase.calculators.orca import ORCA 
from ase.optimize.lbfgs import LBFGS

from src.metal_dock.logger import MetalDockLogger

class ORCAEngine:
    def __init__(self, par, metal_complex):
        self.par = par
        self.metal_complex = metal_complex
        self.logger = MetalDockLogger()

    def run(self):
        """
        Run the ORCA engine.
        """
        if self.par.geom_opt:
            self._run_geom_opt()
        else:
            self._run_single_point()

        return os.getcwd(), self.energy

    def _run_geom_opt(self):
        """
        Run the geometry optimization.
        """
        self.logger.info("STARTING ORCA GEOMETRY OPTIMIZATION...")
        os.makedirs(self.par.output_dir / 'QM' / 'geom_opt', exist_ok=True)

        shutil.copy(self.par.xyz_file, self.par.output_dir / 'QM' / 'geom_opt' / self.par.xyz_file)

        out_file = self.par.output_dir / 'QM' / 'geom_opt' / 'geom.out' 
        if not os.path.exists(out_file):
            self._orca_geom_opt()

        if self._orca_opt_converged(out_file):
            self.energy = self._orca_extract_energy(out_file)
            self.metal_complex.energy = self.energy

            self._orca_extract_charges('geom_opt', out_file)
        
        os.chdir(self.par.output_dir)
        self.metal_complex.create_mol_graph('geom_opt')
        self.metal_complex.add_charges_to_graph()
        self._orca_extract_bond_orders(out_file)

    def _run_single_point(self):
        """
        Run the single point calculation.
        """
        os.makedirs(self.par.output_dir / 'QM' / 'single_point', exist_ok=True)

        shutil.copy(self.par.xyz_file, self.par.output_dir / 'QM' / 'single_point' / self.par.xyz_file)

        out_file = self.par.output_dir / 'QM' / 'single_point' / 'single_point.out'
        if not os.path.exists(out_file):
            self._orca_sp()

        if self._orca_sp_converged(out_file):
            self.energy = self._orca_extract_energy(out_file)
            self.metal_complex.energy = self.energy

            self._orca_extract_charges('single_point', out_file)

        os.chdir(self.par.output_dir)
        self.metal_complex.create_mol_graph('single_point')
        self.metal_complex.add_charges_to_graph()
        self._orca_extract_bond_orders(out_file)

    def _orca_extract_energy(self, log_file: str):
        """
        Extract the energy from the log file.

        Args:
            log_file (str): The path to the log file.

        Returns:
            float: The energy from the log file.
        """
        with open(log_file, 'r') as fin:
            for line in fin:
                if line.startswith('FINAL'):
                    energy = line.split()[4]
                    return energy

    def _orca_extract_charges(self, run_type: str, log_file: str):
        """
        Extract the CM5 charges from the log file.

        Args:
            run_type (str): The type of run (either 'geom_opt' or 'single_point').
            log_file (str): The path to the log file.
        """
        cm5_file = self.par.output_dir / 'QM' / run_type / 'CM5_charges.csv'

        a0, radii, periodic_table = self._orca_load_CM5_model()
        data = self._orca_get_CM5_logfile(log_file, periodic_table, radii)
        qcm5 = self._orca_HirshfeldToCM5(data, a0)

        charges = {}
        for idx, row in qcm5.iterrows():
            charges[f'{row["ATOM"].upper()}{idx+1}'] = row['QCM5']

        # save qm5 to csv
        qcm5.to_csv(cm5_file, index=False, float_format='%6.4f')

        self.metal_complex.charges = charges

    def _orca_HirshfeldToCM5(self, df: pd.DataFrame, a0: pd.DataFrame): 
        """
        Convert the Hirshfeld charges to CM5 charges.

        Args:
            df (pd.DataFrame): DataFrame containing atomic coordinates and atomic numbers.
            a0 (pd.DataFrame): DataFrame containing the A values from the CM5 model.

        Returns:
            pd.DataFrame: DataFrame containing the CM5 charges.
        """
        DVALS = self._orca_get_avals(a0)
        cm5_charges = []
        alpha = 2.474
        for i, r in df.iterrows(): 
            qcm5 = r.QHir
            for j, p in df.iterrows(): 
                if (r.AtNum != p.AtNum): 
                    dist = self._orca_get_distance([r.X, r.Y, r.Z], [p.X, p.Y, p.Z])
                    factor = np.exp(-1.0 * alpha * (dist - r.RAD - p.RAD))
                    qcm5 = qcm5 + factor * DVALS[r.AtNum - 1, p.AtNum - 1]
            cm5_charges.append(qcm5)
        df['QCM5'] = np.array(cm5_charges)

        # Use the new local environment calculation
        radius = 2.0  # Define the radius
        df['LocalEnv'] = [tuple(self._orca_calculate_local_environment(df, atomNum, radius)) for atomNum in df.index]
        uniq_envs = list(set(df['LocalEnv']))
        df['QCM5_AVG'] = [df[df['LocalEnv'] == env].QCM5.mean() for env in df['LocalEnv']]
        df['1.20*CM5'] = df.QCM5_AVG * 1.20
        return df
    
    def _orca_calculate_local_environment(self, df: pd.DataFrame, atom_index: int, radius: float):
        """
        Calculate the local environment of an atom based on atomic numbers within a given radius.

        Args:
            df (pd.DataFrame): DataFrame containing atomic coordinates and atomic numbers.
            atom_index (int): Index of the target atom.
            radius (float): Radius within which to consider neighboring atoms.

        Returns:
            list: List of atomic numbers within the specified radius.
        """
        target_atom = df.iloc[atom_index]
        target_coords = np.array([target_atom['X'], target_atom['Y'], target_atom['Z']])
        local_environment = []

        for i, row in df.iterrows():
            if i == atom_index:
                continue
            neighbor_coords = np.array([row['X'], row['Y'], row['Z']])
            distance = np.linalg.norm(target_coords - neighbor_coords)
            if distance <= radius:
                local_environment.append(row['AtNum'])

        # Optionally sort the atomic numbers for consistency
        local_environment.sort()
        return local_environment
    
    def _orca_get_distance(self, a: list, b: list):
        """
        Get the distance between two points.

        Args:
            a (list): The first point.
            b (list): The second point.

        Returns:
            float: The distance between the two points.
        """
        return(np.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2))

    def _orca_load_CM5_model(self):
        """
        Load the CM5 model.
        """
        cm5_path = os.path.join(os.environ['ROOT_DIR'], 'metal_dock', 'cm5pars.json')
        with open(cm5_path) as tweetfile:
            cm5_model = json.loads(tweetfile.read())
        a0_df = pd.DataFrame.from_dict(cm5_model['A0'])
        rd_df = pd.DataFrame.from_dict(cm5_model['radii'])
        pt_df = pd.DataFrame.from_dict(cm5_model['PeriodicTable'])
        return (a0_df,rd_df,pt_df)
    
    def _orca_get_CM5_logfile(self, log_file: str, pt_df: pd.DataFrame, rd_df: pd.DataFrame):
        """
        Extract the CM5 charges from the log file.

        Args:
            log_file (str): The path to the log file.
            pt_df (pd.DataFrame): DataFrame containing the periodic table.
            rd_df (pd.DataFrame): DataFrame containing the radii.

        Returns:
            pd.DataFrame: DataFrame containing the CM5 charges.
        """
        pt_df["symbol"] = pt_df["symbol"].map(str.strip)
        sym2num = pt_df.set_index(['symbol'])['atomicNumber'].to_dict() 
        num2rad = rd_df.set_index(['RAD_NO'])['VALUE'].to_dict()
        xyz_data = []
        charge_data = []
        data = open(log_file).readlines()
        id_charges = False
        id_coos = False
        for line in data: 
            if 'CARTESIAN COORDINATES (ANGSTROEM)' in line: 
                id_coos = True 
            elif 'CARTESIAN COORDINATES (A.U.)' in line:
                id_coos = False 
            if 'HIRSHFELD ANALYSIS' in line: 
                id_charges=True
            elif 'TIMINGS' in line: 
                id_charges = False
            if id_charges:charge_data.append(line.strip().split()) 
            if id_coos:xyz_data.append(line.strip().split()) 
        hirCharges = pd.DataFrame(charge_data[7:-4],columns=['N','ATOM','QHir','Spin'])
        hirCharges[['N','QHir','Spin']] = hirCharges[['N','QHir','Spin']].apply(pd.to_numeric)
        hirCharges = hirCharges[['N','QHir']]
        xyzcoos = pd.DataFrame(xyz_data[2:-2],columns=['ATOM','X','Y','Z'])
        xyzcoos[['X','Y','Z']] = xyzcoos[['X','Y','Z']].apply(pd.to_numeric)
        final_data = (pd.concat([xyzcoos,hirCharges],axis=1))
        final_data['AtNum'] = [sym2num[s] for s in final_data.ATOM] 
        final_data['RAD'] = [num2rad[s] for s in final_data.AtNum] 
        return(final_data)
    
    def _orca_get_avals(self, a0_df: pd.DataFrame):
        """
        Get the A values from the CM5 model.

        Args:
            a0_df (pd.DataFrame): DataFrame containing the A values from the CM5 model.

        Returns:
            np.ndarray: Array containing the A values.
        """
        num2a = (a0_df.set_index(['A0_NO'])['VALUE'].to_dict())
        list_keys = list(num2a.keys())
        dvals = (np.empty([np.max(list_keys),np.max(list_keys)]))
        for i in range(dvals.shape[0]):
            for j in range(dvals.shape[1]):
                if (i != j) : dvals[i,j] = num2a[i+1]-num2a[j+1]
        dvals[ 0, 5]= 0.0502
        dvals[ 0, 6]= 0.1747
        dvals[ 0, 7]= 0.1671
        dvals[ 5, 6]= 0.0556
        dvals[ 5, 7]= 0.0234
        dvals[ 6, 7]=-0.0346
        # Repitition of the above cooefficients with a negative sign
        dvals[ 5, 0]=-0.0502
        dvals[ 6, 0]=-0.1747
        dvals[ 7, 0]=-0.1671
        dvals[ 6, 5]=-0.0556
        dvals[ 7, 5]=-0.0234
        dvals[ 7, 6]= 0.0346
        return dvals

    def _orca_extract_bond_orders(self, log_file: str):
        """
        Extract the bond orders from the log file.

        Args:
            log_file (str): The path to the log file.
        """
        bond_orders = {}

        with open(log_file, 'r') as file:
            content = file.read()

        pattern = r"Mayer bond orders larger than \d+\.\d+\n(.*?)\n\n"
        match = re.search(pattern, content, re.DOTALL)

        if match:
            bond_order_section = match.group(1)

            bond_pattern = r"B\(\s*(\d+)-\w+\s*,\s*(\d+)-\w+\s*\)\s*:\s*(\d+\.\d+)"
            bonds = re.findall(bond_pattern, bond_order_section)

            for bond in bonds:
                atom1 = int(bond[0])  # Convert to zero-based index
                atom2 = int(bond[1])  # Convert to zero-based index
                order = float(bond[2])

                # Sort atoms
                atoms = tuple(sorted((atom1, atom2)))

                # Store bond order for both directions
                bond_orders[atoms] = order

        # Update the metal_complex graph with bond orders
        for (atom1, atom2), order in bond_orders.items():
            if self.metal_complex.graph.has_edge(atom1, atom2):
                self.metal_complex.graph[atom1][atom2]['bond_order'] = order

        # Remove edges without bond orders from the graph
        for edge in list(self.metal_complex.graph.edges()):
            sorted_edge = tuple(sorted(edge))
            if sorted_edge not in bond_orders:
                self.metal_complex.graph.remove_edge(*edge)

    def _orca_opt_converged(self, log_file: str):
        """
        Check if the geometry optimization has converged.

        Args:
            log_file (str): The path to the log file.

        Returns:
            bool: True if the geometry optimization has converged, False otherwise.
        """
        with open(log_file) as log:
            if 'SUCCESS' in log.read():
                self.logger.info('GEOMETRY CONVERGED')
                return True
            else:
                self.logger.info('GEOMETRY NOT CONVERGED - VERIFY THE PROBLEM IN geom.out AND DELETE FILE BEFORE RUNNING AGAIN')
                sys.exit()

    def _orca_sp_converged(self, log_file):
        """
        Check if the single point calculation has converged.

        Args:
            log_file (str): The path to the log file.

        Returns:
            bool: True if the single point calculation has converged, False otherwise.
        """
        with open(log_file) as log:
            if 'SUCCESS' in log.read():
                self.logger.info('SINGLE POINT CONVERGED')
                return True
            else:
                self.logger.info('SINGLE POINT NOT CONVERGED - VERIFY THE PROBLEM IN single_point.out AND DELETE FILE BEFORE RUNNING AGAIN')
                sys.exit()

    def _orca_geom_opt(self):
        """
        Run the geometry optimization.
        """
        os.chdir(self.par.output_dir / 'QM' / 'geom_opt')
        M = 2 * (self.par.spin * 0.5) + 1

        mol = read(self.par.xyz_file)
        mol.calc = ORCA(label='geom',
                        charge=self.par.charge,
                        mult=M,
                        orcasimpleinput=f'{self.par.orcasimpleinput}',
                        orcablocks=f'%pal nprocs {str(self.par.ncpu)} end %output Print[P_hirshfeld] 1 end {self.par.orcablocks}',
                        )

        opt = LBFGS(mol)
        opt.run(fmax=0.05)

        mol.get_potential_energy()
        mol.write(self.par.output_dir / 'QM' / 'geom_opt' / 'output.xyz')

    def _orca_sp(self):
        """
        Run the single point calculation.
        """
        os.chdir(self.par.output_dir / 'QM' / 'single_point')   
        M = 2 * (self.par.spin * 0.5) + 1

        mol = read(self.par.xyz_file)
        mol.calc = ORCA(label='single_point',
                        charge=self.par.charge,
                        mult=M,
                        orcasimpleinput=f'{self.par.orcasimpleinput}',
                        orcablocks=f'%pal nprocs {str(self.par.ncpu)} end %output Print[P_hirshfeld] 1 end {self.par.orcablocks}',
                        )

        mol.get_potential_energy()
        mol.write(self.par.output_dir / 'QM' / 'single_point' / 'output.xyz')