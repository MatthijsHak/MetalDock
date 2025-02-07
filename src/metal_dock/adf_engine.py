import os
import sys
import re
import subprocess as sp
import scm.plams as scm

from src.metal_dock.logger import MetalDockLogger

class ADFEngine:
    def __init__(self, par, metal_complex):
        self.par = par
        self.metal_complex = metal_complex
        self.logger = MetalDockLogger()

    def run(self):
        """
        Run the ADF engine.
        """
        if self.par.geom_opt:
            return self._run_geom_opt()
        else:
            return self._run_single_point()

    def _run_geom_opt(self):
        """
        Run the geometry optimization.
        """
        os.makedirs(self.par.output_dir / 'QM', exist_ok=True)

        log_file = self.par.output_dir / 'QM' / 'geom_opt' / 'plamsjob' / 'ams.log'
        if not os.path.isdir(self.par.output_dir / 'QM' / 'geom_opt' / 'plamsjob'):
            self._adf_geom_opt()

        if self._adf_opt_converged(log_file):
            energy = self._adf_extract_energy(log_file)
            self.metal_complex.energy = energy

            self._adf_extract_CM5('geom_opt')
            self._adf_extract_charges('geom_opt')

        os.chdir(self.par.output_dir)
        self.metal_complex.create_mol_graph('geom_opt')
        self.metal_complex.add_charges_to_graph()
        self._adf_extract_bond_orders('geom_opt')

    def _run_single_point(self):
        """
        Run the single point calculation.
        """
        os.makedirs(self.par.output_dir/ 'QM', exist_ok=True)

        log_file = self.par.output_dir / 'QM' / 'single_point' / 'plamsjob' / 'ams.log'
        if not os.path.isdir(self.par.output_dir/ 'QM' / 'single_point' / 'plamsjob'):
            self._adf_sp()

        if self._adf_sp_converged(log_file):
            energy = self._adf_extract_energy(log_file)
            self.metal_complex.energy = energy

            self._adf_extract_CM5('single_point')
            self._adf_extract_charges('single_point')

        os.chdir(self.par.output_dir)
        self.metal_complex.create_mol_graph('single_point')
        self.metal_complex.add_charges_to_graph()
        self._adf_extract_bond_orders('single_point')

    def _adf_extract_energy(self, log_file):
        """
        Extract the energy from the log file.
        """
        with open(log_file, 'r') as fin:
            for line in fin:
                if 'kcal/mol' in line:
                    energy = line.split()[4]
                    return energy
                
    def _adf_extract_CM5(self, run_type):
        """
        Extract the CM5 charges from the log file.
        """
        rkf_file = self.par.output_dir / 'QM' / run_type / 'plamsjob' / 'adf.rkf'
        charges_file = self.par.output_dir / 'QM' / run_type / 'plamsjob' / 'CM5_charges'
        amsreport_path = os.path.join(os.environ['AMSBIN'], 'amsreport')
        sp.call([amsreport_path + ' ' + str(rkf_file) + ' CM5 > ' + str(charges_file)], shell=True)

    def _adf_extract_charges(self, run_type):
        """
        Extract the CM5 charges from the log file.
        """
        charges = {}
        with open(self.par.output_dir / 'QM' / run_type / 'plamsjob' / 'CM5_charges', 'r') as file:
            for idx, line in enumerate(file):
                # skip first line
                if idx == 0:
                    continue
                line = line.strip().split()
                # if empty line skip
                if len(line) == 0:
                    continue
                charges[line[0]] = line[-1]
        charges = {k.replace('(', '').replace(')', '').upper(): v for k, v in charges.items()}
        self.metal_complex.charges = charges

    def _adf_extract_bond_orders(self, run_type):
        """
        Extract the bond orders from the log file.
        """
        log_file = self.par.output_dir / 'QM' / run_type / 'plamsjob' / 'plamsjob.out'
        
        bond_orders = {}
        with open(log_file, 'r') as file:
            content = file.read()
        
        # Define the pattern to capture the bond order section
        pattern = r"Description: Mayer bond orders\nOnly bonds with bond orders > 0\.200 are printed\.\n\n Index  Atom    Index  Atom    BondOrder\n(.*?)\n\n"
        match = re.search(pattern, content, re.DOTALL)
        
        if match:
            bond_order_section = match.group(1)
            
            # Define the pattern to extract bond order information
            bond_pattern = r"\s*(\d+)\s+\w+\s+(\d+)\s+\w+\s+([\d.]+)"
            bonds = re.findall(bond_pattern, bond_order_section)
            
            for bond in bonds:
                atom1 = int(bond[0])-1
                atom2 = int(bond[1])-1
                order = float(bond[2])

                # Sort atoms
                atoms = tuple(sorted((atom1, atom2)))
                # Store bond order for both directions
                bond_orders[atoms] = order


        # Update the metal_complex graph with bond orders
        for (atom1, atom2), order in bond_orders.items():
            if self.metal_complex.graph.has_edge(atom1, atom2):
                self.metal_complex.graph[atom1][atom2]['bond_order'] = order

        # if an edge does not have a bond order then we need to remove it from the graph
        for edge in self.metal_complex.graph.edges():
            # Sort the edge to ensure order independence
            sorted_edge = tuple(sorted(edge))
            if sorted_edge not in bond_orders:
                self.metal_complex.graph.remove_edge(*edge)

    def _adf_opt_converged(self, ams_log):
        """
        Check if the geometry optimization has converged.
        """
        with open(ams_log) as log:
            if 'Geometry optimization converged' in log.read():
                self.logger.info('GEOMETRY CONVERGED')
                return True
            else:
                self.logger.info('GEOMETRY NOT CONVERGED\nPLEASE CHECK THE ams.log FILE IN THE QM DIRECTORY')
                return sys.exit()

    def _adf_sp_converged(self, ams_log):
        """
        Check if the single point calculation has converged.
        """
        with open(ams_log) as log:
            if 'NORMAL TERMINATION' in log.read():
                self.logger.info('SINGLE POINT CONVERGED\n')
                return True
            else:
                self.logger.info('SINGLE POINT NOT CONVERGED\nPLEASE CHECK THE ams.log FILE IN THE QM DIRECTORY')
                return sys.exit()

    def _adf_geom_opt(self):
        """
        Run the geometry optimization.
        """
        scm.init(folder=self.par.output_dir / 'QM' / 'geom_opt')
        m = scm.Molecule(self.par.xyz_file)
        m.properties.charge = str(self.par.charge)

        s = scm.Settings()
        s.input.ams.Task = 'GeometryOptimization'
        s.input.ams.properties.bondorders = 'Yes'
        s.input.adf.bondorders.TypeForAMS = 'Mayer'
        s.input.adf.scf.iterations = '500'
        s.input.adf.AtomicChargesTypeForAMS = 'CM5'
        s.input.adf.basis.type = self.par.basis_set.upper()
        s.input.adf.basis.core = 'None'

        self._set_functional(s)
        self._set_dispersion(s)
        self._set_spin(s)
        self._set_relativity(s)
        self._set_solvent(s)

        j = scm.AMSJob(molecule=m, settings=s)
        result = j.run()
        j.results.get_main_molecule().write(self.par.output_dir / 'QM' / 'geom_opt' / 'output.xyz', 'xyz')
        scm.finish()

    def _adf_sp(self):
        """
        Run the single point calculation.
        """
        scm.init(folder=self.par.output_dir / 'QM' / 'single_point')
        m = scm.Molecule(self.par.xyz_file)
        m.properties.charge = str(self.par.charge)

        s = scm.Settings()
        s.input.ams.Task = 'SinglePoint'
        s.input.ams.properties.bondorders = 'Yes'
        s.input.adf.bondorders.TypeForAMS = 'Mayer'
        s.input.adf.scf.iterations = '500'
        s.input.adf.AtomicChargesTypeForAMS = 'CM5'
        s.input.adf.basis.type = self.par.basis_set.upper()
        s.input.adf.basis.core = 'None'

        self._set_functional(s)
        self._set_dispersion(s)
        self._set_spin(s)
        self._set_relativity(s)
        self._set_solvent(s)

        j = scm.AMSJob(molecule=m, settings=s)
        result = j.run()
        scm.finish()
        j.results.get_main_molecule().write(self.par.output_dir / 'QM' / 'single_point' / 'output.xyz', 'xyz')
        scm.finish()

    def _set_functional(self, s):
        """
        Set the functional type.
        """
        if self.par.functional_type.lower() == '':
            raise ValueError('Functional type not specified. Please specify functional type (LDA, GGA, METAGGA, HYBRID, METAHYBRID) in the input file')

        if self.par.functional_type.lower() == 'lda':
            s.input.adf.xc.lda = self.par.functional.upper()

        if self.par.functional_type.lower() == 'metagga':
            s.input.adf.xc.metagga = self.par.functional.upper()

        if self.par.functional_type.lower() == 'gga':
            s.input.adf.xc.gga = self.par.functional.upper()

        if self.par.functional_type.lower() == 'hybrid':
            s.input.adf.xc.hybrid = self.par.functional.upper()

        if self.par.functional_type.lower() == 'metahybrid':
            s.input.adf.xc.metahybrid = self.par.functional.upper()

    def _set_dispersion(self, s):
        """
        Set the dispersion type.
        """
        if self.par.dispersion is not None:
            s.input.adf.xc.dispersion = self.par.dispersion.upper()

    def _set_spin(self, s):
        """
        Set the spin type.
        """
        if self.par.spin != 0:
            s.input.adf.unrestricted = 'yes'
            s.input.adf.spinpolarization = str(self.par.spin)

    def _set_relativity(self, s):
        """
        Set the relativity type.
        """
        if self.par.relativity == 'ZORA':
            s.input.adf.relativity.formalism = 'ZORA'
            s.input.adf.relativity.level = 'Scalar'

    def _set_solvent(self, s):
        """
        Set the solvent type.
        """
        if self.par.solvent != '':
            s.input.adf.Solvation.Solv = f"Name={self.par.solvent}"