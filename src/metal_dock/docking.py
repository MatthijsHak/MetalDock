import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

from scipy.spatial.distance import cdist
from src.metal_dock.logger import MetalDockLogger

class Docking:
    def __init__(self, par, metal_complex, protein):
        self.logger = MetalDockLogger() 
        self.par = par
        self.metal_complex = metal_complex
        self.protein = protein
        self.box_centre, self.box_size = self._get_box_parameters()

    def run(self):
        """
        Run the MetalDock process.
        """
        self.logger.info("STARTING METALDOCK DOCKING PROCESS...")
        # move the ligand and receptor pdbqt files to the docking directory
        shutil.move(self.par.output_dir / 'file_prep' / f'{self.par.name_ligand}.pdbqt', 
                   self.par.output_dir / 'docking' / f'{self.par.name_ligand}.pdbqt')
        shutil.move(self.par.output_dir / 'file_prep' / f'clean_{self.par.name_protein}.pdbqt', 
                   self.par.output_dir / 'docking' / f'{self.par.name_protein}.pdbqt')

        gpf_path = self.par.output_dir / 'docking' / f'{self.par.name_ligand}_{self.par.name_protein}.gpf'
        dpf_path = self.par.output_dir / 'docking' / f'{self.par.name_ligand}_{self.par.name_protein}.dpf'

        if self.par.parameter_file == 'metal_dock.dat':
            parameter_path = os.environ['ROOT_DIR']+'/metal_dock/'+self.par.parameter_file
        else:
            parameter_path = self.par.parameter_file

        ligand_pdbqt_path = self.par.output_dir / 'docking' / f'{self.par.name_ligand}.pdbqt'
        receptor_pdbqt_path = self.par.output_dir / 'docking' / f'{self.par.name_protein}.pdbqt'

        if not gpf_path.exists():
            self._create_gpf_file(ligand_pdbqt_path, receptor_pdbqt_path, gpf_path, parameter_path)
        else:
            os.remove(gpf_path)
            self._create_gpf_file(ligand_pdbqt_path, receptor_pdbqt_path, gpf_path, parameter_path)

        if not dpf_path.exists():
            self._create_dpf_file(dpf_path, gpf_path, parameter_path)
        else:
            os.remove(dpf_path)
            self._create_dpf_file(dpf_path, gpf_path, parameter_path)

        dlg_path = self.par.output_dir / 'docking' / f'{self.par.name_ligand}_{self.par.name_protein}.dlg'
        dock_dir_path = self.par.output_dir / 'docking'
        
        results_dir_path = self.par.output_dir / 'results'
        results_dir_path.mkdir(exist_ok=True)
        autogrid_logfile = self.par.output_dir / 'docking' / 'autogrid.log'
        autodock_logfile = self.par.output_dir / 'docking' / 'autodock.log'

        os.chdir(dock_dir_path)
        self._run_autogrid(autogrid_logfile, gpf_path)
        self._run_autodock(autodock_logfile, dpf_path)
        if not self.par.score_only:
            self._write_conformations(dlg_path, dock_dir_path)
            # self._adding_and_optimizing_hydrogens()
            self._clean_dummy_atoms()
            self._write_pdbqt_to_xyz()
            self._write_pose_to_pdb()

    def _write_pdbqt_to_xyz(self):
        """
        Write the pdbqt file to an xyz file
        """
        docking_dir = self.par.output_dir / 'docking'
        results_dir = self.par.output_dir / 'results'

        for i in range(1, self.par.num_poses+1):
            pose_dir = results_dir / f'pose_{i}'
            pose_dir.mkdir(exist_ok=True)

            pdqt_in = docking_dir / f'{self.par.name_ligand}_{i}.pdbqt'
            pdqt_out = pose_dir / f'{self.par.name_ligand}_{i}.xyz'
            self._write_pose_to_xyz(pdqt_out, pdqt_in)

    def _write_pose_to_xyz(self, xyz_file, pdbqt_file):
        """
        Write the pose to an xyz file
        """
        output_lines = []
        with open(pdbqt_file, 'r') as fin:
            for line in fin:
                if 'ATOM' in line or 'HETATM' in line:
                    splits = line.strip().split()
                    output_lines.append(f'{splits[2]:>2} {float(splits[6]):>8.3f} {float(splits[7]):>8.3f} {float(splits[8]):>8.3f}\n')

        # insert the total number of atoms at the top of the file
        output_lines.insert(0, f'{len(output_lines)}\n\n')

        with open(xyz_file, 'w') as fout:
            for line in output_lines:
                fout.write(line)

    def analyze_results(self):
        # copy the dlg file to results directory
        dlg_in = self.par.output_dir / 'docking' / f'{self.par.name_ligand}_{self.par.name_protein}.dlg'
        dlg_out = self.par.output_dir / 'results' / 'docking_results.dlg'
        shutil.copyfile(dlg_in, dlg_out)

        # copy the pdb file of the receptor to results directory
        pdb_in = self.par.output_dir / 'file_prep' / f'clean_{self.par.name_protein}.pdb'
        pdb_out = self.par.output_dir / 'results' / f'clean_{self.par.name_protein}.pdb'
        shutil.copyfile(pdb_in, pdb_out)

        binding_energies, binding_efficiencies = self._extract_dlg(dlg_out)
        self.logger.info('\n#==============================================================================#')
        self.logger.info('DOCKING RESULTS:')
        self.logger.info(f'#==============================================================================#\n')
        self.logger.info('Ligand Efficiency = (binding energy) / (number of heavy atoms in metal complex)')
        self.logger.info('Interacting Residues = residues within 4 Angstrom of the metal complex\n')

        if self.par.score_only:
            self.logger.info(f'Binding Energy: {binding_energies[0]:7.4f} kcal/mol')
            self.logger.info('\n')
        else:
            for i in range(self.par.num_poses):
                self.logger.info(f"Pose {i+1}:")
                self.logger.info("-------------")
                pose_path = self.par.output_dir / 'results' / f'pose_{i+1}' / f'{self.par.name_ligand}_{i+1}.xyz'
                pdb_path = self.par.output_dir / 'results' / f'clean_{self.par.name_protein}.pdb'
                pose_residues = self._extract_interacting_residues(pose_path, pdb_path)
                self.logger.info(f'Binding Energy: {binding_energies[i]:7.4f} kcal/mol')
                self.logger.info(f'Ligand Efficiency: {binding_efficiencies[i]:7.4f} kcal/mol')
                self.logger.info(f'Interacting Residues:')
                for residue in pose_residues:
                    self.logger.info(f'Residue: {residue[0]}, ID: {residue[1]:>3}')
                self.logger.info('\n')

        if self.par.rmsd:
            rmsd_path = self.par.output_dir / 'file_prep' / f'{self.par.name_ligand}_c.xyz'
            results_dir_path = self.par.output_dir / 'results'
            print_list, rmsd_list, avg_rmsd, stdv_rmsd, var_rmsd = self._calculate_rmsd(rmsd_path, results_dir_path)

            for line in print_list:
                self.logger.info(line)
        self.logger.info("THE PRINTED POSES AND PROTEIN CAN BE FOUND IN THE RESULTS DIRECTORY")
        self.logger.info("EACH PDB FILE IN THE RESULTS/POSE_X DIRECTORY CAN BE VISUALIZED WITH E.G. PYMOL")

        self.logger.info('\n#==============================================================================#')
        self.logger.info("METALDOCK SUCCESFULLY COMPLETED")
        self.logger.info('#==============================================================================#\n')

    def _extract_interacting_residues(self, pose, protein, cutoff=4.0):
        pose_xyz = []

        with open(pose, 'r') as fin:
            for _ in range(2):
                next(fin)
            for line in fin:
                split = line.strip().split()[1:4]
                pose_xyz.append([float(split[0]), float(split[1]), float(split[2])])

        interacting_residues = []
        interacting_residues_xyz = []
        # calculate all the residues that are within 4 angstroms of the ligand
        with open(protein,'r') as fin:
            for line in fin:
                if 'ATOM' in line or 'HETATM' in line:
                    # calculate the distance between the pose_xyz and the 
                    split = line.strip().split()
                    interacting_residues.append((split[3], split[5]))
                    interacting_residues_xyz.append([float(split[6]), float(split[7]), float(split[8])])

        distance_matrix = cdist(pose_xyz, interacting_residues_xyz,'euclidean')

        # Find indices of residues within 4 angstroms
        within_4_angstrom_indices = np.any(distance_matrix <= 4.0, axis=0)

        # Get the unique residues that are within 4 angstroms
        unique_residues_within_4_angstrom = set()
        for i, within in enumerate(within_4_angstrom_indices):
            if within:
                unique_residues_within_4_angstrom.add(interacting_residues[i])

        # Convert set to a list for further processing or output
        unique_residues_within_4_angstrom = list(unique_residues_within_4_angstrom)

        # order based on the id of the residues 
        unique_residues_within_4_angstrom = sorted(unique_residues_within_4_angstrom, key=lambda x: int(x[1]))

        return unique_residues_within_4_angstrom

    def _extract_dlg(self, dlg_file):
        """Extracts all occurrences of the Estimated Free Energy of Binding from a docking output file."""
        binding_energies = []
        
        with open(dlg_file, 'r') as file:
            content = file.read()  # Read the entire file as a string

        # Find all occurrences of "Estimated Free Energy of Binding" followed by a number
        matches = re.findall(r'Estimated Free Energy of Binding\s*=\s*([-\d.]+)\s*kcal/mol', content)

        # Convert matches to float values
        binding_energies = [float(value) for value in matches]

        # only take the first n poses 
        binding_energies = binding_energies[:self.par.num_poses]
        binding_efficiencies = [binding_energy / self.par.n_heavy_atoms for binding_energy in binding_energies]
        return binding_energies, binding_efficiencies

    def _write_pose_to_pdb(self, residue_name='UNK'):
        """
        Write the pose to a pdb file

        Args:
            residue_name (str): The three letter code for thet metal complex residue
        """
        results_dir = self.par.output_dir / 'results'
        for i in range(1, self.par.num_poses+1):
            pose_dir = results_dir / f'pose_{i}'
            xyz_file = pose_dir / f'{self.par.name_ligand}_{i}.xyz'
            pdb_file = pose_dir / f'{self.par.name_ligand}_{i}.pdb'
            # read the xyz file 
            atoms = []
            with open(xyz_file, 'r') as fin:
                for _ in range(2):
                    next(fin)
                for line in fin:
                    split = line.strip().split()
                    atoms.append([split[0], [float(split[1]), float(split[2]), float(split[3])]])

            # atom index
    def _write_pose_to_pdb(self, residue_name='UNK'):
        """
        Write the pose to a pdb file

        Args:
            pdb_file (str): The path to the PDB file to write the pose to.
        """
        results_dir = self.par.output_dir / 'results'
        for i in range(1, self.par.num_poses+1):
            pose_dir = results_dir / f'pose_{i}'
            xyz_file = pose_dir / f'{self.par.name_ligand}_{i}.xyz'
            pdb_file = pose_dir / f'{self.par.name_ligand}_{i}.pdb'
            # read the xyz file 
            atoms = []
            with open(xyz_file, 'r') as fin:
                for _ in range(2):
                    next(fin)
                for line in fin:
                    split = line.strip().split()
                    atoms.append([split[0], [float(split[1]), float(split[2]), float(split[3])]])

            # atom index
            atom_index = 1
            with open(pdb_file, 'w') as f:
                for atom in atoms:
                    atom_type = f'{atom[0]}{atom_index-1}'
                    f.write(f"HETATM{atom_index:>5} {atom_type:>3}  {residue_name} A   1    {atom[1][0]:>8.3f}{atom[1][1]:>8.3f}{atom[1][2]:>8.3f}  1.00  0.00          {atom[0]:>2}\n")
                    atom_index += 1

                for edge in self.metal_complex.graph.edges():
                        f.write(f"CONECT {self.metal_complex.atom_index_mapping[edge[0]]['pdbqt_index']+1:>4} {self.metal_complex.atom_index_mapping[edge[1]]['pdbqt_index']+1:>4}\n")
                f.write('ENDMDL\n')

    def _write_xyz_file(self, pdbqt_lines, xyz_file):
        with open(xyz_file, 'w') as fout:
            fout.write(f'{len(pdbqt_lines)}\n')
            fout.write(f'{self.par.name_ligand}\n')
            for line in pdbqt_lines:
                fout.write(f'{line[1]:<2} {float(line[2]):>8.3f} {float(line[3]):>8.3f} {float(line[4]):>8.3f}\n')

    def _clean_dummy_atoms(self):
        """
        Remove dummy atoms from the metal complex graph and the pdbqt file.
        """
        # Collect nodes to remove in a separate list
        nodes_to_remove = [node for node in self.metal_complex.graph.nodes() if self.metal_complex.graph.nodes[node]['element'] == 'DD']

        # Remove nodes after iteration
        for node in nodes_to_remove:
            self.metal_complex.graph.remove_node(node)

        # remove dummy atoms from the pdbqt file
        for n in range(self.par.num_poses):
            self._delete_dummy_atom(self.par.output_dir / 'docking' / f'{self.par.name_ligand}_{n+1}.pdbqt')

    def _delete_dummy_atom(self, pdbqt_file):
        """
        Delete dummy atoms from the pdbqt file.

        Args:
            pdbqt_file (str): The path to the PDBQT file to delete dummy atoms from.
        """
        with open(pdbqt_file,'r') as fin:
            with open('output.pdbqt','w') as fout:
                for line in fin:
                    if 'DD' in line:
                        pass
                    else:
                        fout.write(line)
        shutil.move('output.pdbqt', pdbqt_file)

    def _write_conformations(self, dlg_path, output_path):
        """
        Write the pose to a pdbqt file
        """
        with open(dlg_path, 'r') as fin:
            pose = []
            atom_id = 0
            mol_id = 1
            in_docked_block = False  
            docked_block = []  
            for line in fin:
                if 'DOCKED: ROOT' in line:
                    in_docked_block = True
                    docked_block.append(line)
                elif in_docked_block and 'TER' in line:
                    in_docked_block = False
                    docked_block.append(line)
                    with open(output_path / f'{self.par.name_ligand}_{mol_id}.pdbqt', 'w') as output_file:
                        for block_line in docked_block:
                            cleaned_line = block_line.replace('DOCKED: ', '', 1)
                            output_file.write(cleaned_line)
                    mol_id += 1
                    docked_block = []  
                elif in_docked_block:
                    docked_block.append(line)

    def _run_autodock(self, log_file, dpf_path):
        """
        Run the autodock4 program.

        Args:
            dpf_path (str): The path to the DPF file to run.
        """
        autodock4 = os.path.join(os.environ['ROOT_DIR'], 'external', 'AutoDock', 'autodock4')

        with open(log_file, 'w') as log_file:
            subprocess.call(
                [f'{autodock4} -p {dpf_path}'],
                shell=True,
                stdout=log_file,
                stderr=subprocess.STDOUT
            )

    def _run_autogrid(self, log_file, gpf_path):
        """
        Run the autogrid4 program.

        Args:
            gpf_path (str): The path to the GPF file to run.
        """
        autogrid4 = os.path.join(os.environ['ROOT_DIR'],'external','AutoDock','autogrid4')
        with open(log_file, 'w') as log_file:
            subprocess.call([f'{autogrid4} -p {gpf_path}'], 
                            shell=True, 
                            stdout=log_file, 
                            stderr=subprocess.STDOUT)

    def _create_dpf_file(self, dpf_path, gpf_path, parameter_path):
        """
        Create the DPF file for docking.

        Args:
            dpf_path (str): The path to the DPF file to create.
            gpf_path (str): The path to the GPF file to read.
            parameter_path (str): The path to the parameter file to read.
        """
        gpf_file = open(gpf_path,'r')
        gpf_lines = [line.split() for line in gpf_file]

        ligand_type = gpf_lines[5]
        del ligand_type[0]
        del ligand_type[-4:]
        ligand_type_str = ' '.join(ligand_type)

        dpf_file = open(dpf_path,'w')
        dpf_file.write('autodock_parameter_version 4.2       # used by autodock to validate parameter set\n')
        dpf_file.write('parameter_file '+str(parameter_path)+' # parameter library filename\n')
        dpf_file.write('outlev 1                             # diagnostic output level\n')
        dpf_file.write('intelec                              # calculate internal electrostatics\n')
        dpf_file.write('seed pid time                        # seeds for random generator\n')
        dpf_file.write('ligand_types '+ligand_type_str+'             # atoms types in ligand\n')
        dpf_file.write('fld '+self.par.name_protein+'.maps.fld              # grid_data_file\n')
        for i in range(0,len(ligand_type)):
            dpf_file.write('map '+self.par.name_protein+'.'+ligand_type[i]+'.map                 # atom-specific affinity map\n')

        dpf_file.write('elecmap '+self.par.name_protein+'.e.map             # electrostatics map\n')
        dpf_file.write('desolvmap '+self.par.name_protein+'.d.map           # desolvation map\n\n')
        dpf_file.write('move '+self.par.name_ligand+'.pdbqt                # small molecule\n')

        if self.par.random_pos == True:
            dpf_file.write('tran0 random                         # initial coordinates/A or random\n')
            dpf_file.write('quaternion0 random                   # initial orientation\n')
            dpf_file.write('dihe0 random                         # initial dihedrals (relative) or random\n')

        if self.par.ga_dock == True and self.par.sa_dock == False:
            dpf_file.write('# GA parameters\n')
            dpf_file.write('ga_pop_size '+str(self.par.dock_algorithm[0])+'                      # number of individuals in population\n')
            dpf_file.write('ga_num_evals '+str(self.par.dock_algorithm[1])+'                 # maximum number of energy evaluations\n')
            dpf_file.write('ga_num_generations '+str(self.par.dock_algorithm[2])+'             # maximum number of generations\n')
            dpf_file.write('ga_elitism '+str(self.par.dock_algorithm[3])+'                         # number of top individuals to survive to next generation\n')
            dpf_file.write('ga_mutation_rate '+str(self.par.dock_algorithm[4])+'                 # rate of gene mutation\n')
            dpf_file.write('ga_crossover_rate '+str(self.par.dock_algorithm[5])+'                # rate of crossover\n')
            dpf_file.write('ga_window_size '+str(self.par.dock_algorithm[6])+'                    # number of preceding generation when deciding threshold for worst individual current population\n')
            dpf_file.write('ga_cauchy_alpha 0.0                  # Alpha parameter of Cauchy distribution\n')
            dpf_file.write('ga_cauchy_beta 1.0                   # Beta parameter Cauchy distribution\n')

            dpf_file.write('# Local Search Parameters\n')
            dpf_file.write('sw_max_its 300                       # iterations of Solis & Wets local search\n')
            dpf_file.write('sw_max_succ 4                        # consecutive successes before changing rho\n')
            dpf_file.write('sw_max_fail 4                        # consecutive failures before changing rho\n')
            dpf_file.write('sw_rho 1.0                           # size of local search space to sample\n')
            dpf_file.write('sw_lb_rho 0.01                       # lower bound on rho\n')
            dpf_file.write('ls_search_freq 0.06                  # probability of performing local search on individual\n')
            # dpf_file.write('do_local_only 20\n')
            dpf_file.write('# Activate LGA\n')
            dpf_file.write('set_ga                               # set the above parameters for GA or LGA\n')
            dpf_file.write('set_psw1                             # set the above pseudo-Solis & Wets parameters\n')
            dpf_file.write('ga_run '+str(self.par.num_poses)+'                             # do this many hybrid GA-LS runs\n')
        if self.par.ga_dock == False and self.par.sa_dock == True:
            dpf_file.write('# SA Parameters\n')
            dpf_file.write('tstep 2.0\n')
            #dpf_file.write('e0max 0.0 10000                      # max initial energy; max number of retries\n')
            dpf_file.write('linear_schedule                      # linear_schedule or geometric_schedule\n')
            dpf_file.write('rt0 500                              # initial annealing temperature (absolute tmperature multiplied by gas constant\n')
            dpf_file.write('rtrf '+str(self.par.dock_algorithm[0])+'           # annealing temperature reductin factor < 1 cools > 1 heats system\n')
            dpf_file.write('runs '+str(self.par.dock_algorithm[1])+'           # number of docking runs\n')
            dpf_file.write('cycles '+str(self.par.dock_algorithm[2])+'         # number of temperature reduction cycles\n')
            dpf_file.write('accs 30000                           # maximum number of accepted steps per cycle\n')
            dpf_file.write('rejs 30000                           # maximum number of rejected steps per cycle\n')
            dpf_file.write('select m                             # m selects the minimum state, 1 selects the last state during each cycle\n')
            dpf_file.write('trnrf 1.0                            # per cycle reduction factor for translation steps\n')
            dpf_file.write('quarf 1.0                            # per cycle reduction factor for orientation steps\n')
            dpf_file.write('dihrf 1.0                            # per cycle reduction factor for torsional dihedral steps\n')

            dpf_file.write('# Activate SA\n')
            dpf_file.write('simanneal '+str(self.par.num_poses)+'                         # run this many SA docking\n')
        
        elif self.par.score_only == True:
            dpf_file.write('epdb                                # used to calculate the energy of a particular ligand conformation')

        dpf_file.write('analysis                             # perforem a ranked cluster analysis\n')

    def _create_gpf_file(self, ligand_pdbqt_path, receptor_pdbqt_path, gpf_path, parameter_path):
        """
        Create the GPF file for docking.

        Args:
            gpf_path (str): The path to the GPF file to create.
            parameter_path (str): The path to the parameter file to read.
        """
        prepare_gpf4 = os.path.join(os.environ['MGLTOOLS'], 'prepare_gpf4.py')
        command = os.environ['PYTHON_3']+f" {prepare_gpf4} -l {ligand_pdbqt_path} -r {receptor_pdbqt_path} -p parameter_file={parameter_path} -p npts='{self.box_size[0]},{self.box_size[1]},{self.box_size[2]}' -p gridcenter='{self.box_centre[0]:.6},{self.box_centre[1]:.6},{self.box_centre[2]:.6}' -o {gpf_path}"
        subprocess.call([command], shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

        gpf = open(gpf_path,'a')
        gpf.write(f'nbp_r_eps 0.25 23.2135   12 6  NA TZ\n')
        gpf.write(f'nbp_r_eps 2.10  3.8453   12 6  OA Zn\n')
        gpf.write(f'nbp_r_eps 2.25  7.5914   12 6  SA Zn\n')
        gpf.write(f'nbp_r_eps 1.00  0.0000   12 6  HD Zn\n')
        gpf.write(f'nbp_r_eps 2.00  0.0060   12 6  NA Zn\n')
        gpf.write(f'nbp_r_eps 2.00  0.2966   12 6  N  Zn\n')

        if self.par.internal_param == False:
            gpf.write(f'nbp_r_eps 2.20  {self.par.parameter_set[0]:>.4f}   12 10 NA {self.par.metal_symbol}\n')
            gpf.write(f'nbp_r_eps 2.25  {self.par.parameter_set[1]:>.4f}   12 10 OA {self.par.metal_symbol}\n')
            gpf.write(f'nbp_r_eps 2.30  {self.par.parameter_set[2]:>.4f}   12 10 SA {self.par.metal_symbol}\n')
            gpf.write(f'nbp_r_eps 1.00  {self.par.parameter_set[3]:>.4f}   12 6  HD {self.par.metal_symbol}\n')
        gpf.close()

    def _get_box_parameters(self):
        """
        Get the box parameters for docking.
        """
        if self.par.dock_x is not None and self.par.dock_y is not None and self.par.dock_z is not None:
            box_centre = [self.par.dock_x, self.par.dock_y, self.par.dock_z]
        else:
            box_centre = self._get_box_centre()

        # if one value in the list of box size is not 0 then use that value
        if any(x != 0 for x in self.par.box_size) and self.par.scale_factor == 0:
            npts = [x * 2.66 for x in self.par.box_size] # Convert Å to grid point
            if [int(x) for x in npts] == npts:
                box_size =  [int(x) for x in npts]
            else:
                # box_size = math.ceil(npts)
                box_size = [math.ceil(x) for x in npts]
                if self.par.method != 'mc':  # Check if method is not 'mc'
                    self.logger.info('SPACING BETWEEN GRID POINTS IS STANDARD SET TO 0.375 Å')
                    self.logger.info('BOX SIZE MUST BE INTEGER GRID POINTS WHICH WAS NOT FOUND')
                    self.logger.info(f'BOX SIZE SIDE ROUNDED UP AND SET TO {box_size[0] / 2.66:.3f} {box_size[1] / 2.66:.3f} {box_size[2] / 2.66:.3f} Å\n')

        # if all values in the list of box size are 0 then calculate the box size
        if all(x == 0 for x in self.par.box_size) and self.par.scale_factor != 0:
            box_size = self._calculate_box_size(self.par.metal_symbol, 0.375, self.par.scale_factor)

        if self.par.box_size != 0 and self.par.scale_factor != 0:
            if self.par.method != 'mc':  # Check if method is not 'mc'
                self.logger.info("CANNOT SELECT BOXSIZE AND SCALE FACTOR - SET ONE VALUE TO 0")
            sys.exit()

        if self.par.box_size == 0 and self.par.scale_factor == 0:
            if self.par.method != 'mc':  # Check if method is not 'mc'
                self.logger.info("CANNOT SELECT BOXSIZE AND SCALE FACTOR - SET ONE VALUE GREATER THAN 0")
            sys.exit()
        
        return box_centre, box_size
    
    def _get_box_centre(self):
        """
        Get the box centre for docking.

        Args:
            metal_symbol (str): The symbol of the metal to get the centre for.

        Returns:
            list: The coordinates of the box centre.
        """
        dock_x = None
        dock_y = None
        dock_z = None

        # Iterate over all nodes in the graph
        if self.par.method.lower() == 'mc':
            # open the xyz file and read the coordinates
            with open(self.par.output_dir / 'docking' / f'{self.par.name_ligand}_{self.par.name_protein}.xyz', 'r') as f:
                dock_x, dock_y, dock_z = f.read().split()
        else:
            for node, data in self.metal_complex.graph.nodes(data=True):
                if data['element'] == self.par.metal_symbol:
                    dock_x, dock_y, dock_z = data['xyz']
                    break

        if dock_x is None:
            raise ValueError('metal symbol not found in the graph')

        dock = [dock_x, dock_y, dock_z]
        return dock
    
    def _calculate_box_size(self, metal_symbol, spacing, scale_factor):
        """
        Calculate the box size for docking.

        Args:
            metal_symbol (str): The symbol of the metal to calculate the box size for.
            spacing (float): The spacing between grid points.
            scale_factor (float): The scale factor for the box size.

        Returns:
            list: The box size.
        """
        coordinates = []
        x_axis = []
        y_axis = []
        z_axis = []
        metal = None

        # Iterate over all nodes in the graph
        for node, data in self.metal_complex.graph.nodes(data=True):
            element = data['element']
            xyz = data['xyz']
            coordinates.append(xyz)
            x_axis.append(xyz[0])
            y_axis.append(xyz[1])
            z_axis.append(xyz[2])

            if element == metal_symbol:
                metal = np.array(xyz)

        if metal is None:
            raise ValueError('metal symbol not found in the graph')

        # Shift axis to centre at metal
        x_dist = np.abs(np.max(x_axis - metal[0]) - np.min(x_axis - metal[0])) * scale_factor
        y_dist = np.abs(np.max(y_axis - metal[1]) - np.min(y_axis - metal[1])) * scale_factor
        z_dist = np.abs(np.max(z_axis - metal[2]) - np.min(z_axis - metal[2])) * scale_factor

        # Limit distances to a maximum of 20
        x_dist = min(x_dist, 20)
        y_dist = min(y_dist, 20)
        z_dist = min(z_dist, 20)

        # Calculate the number of points
        x_npts = (round(x_dist / spacing)) & (-2)
        y_npts = (round(y_dist / spacing)) & (-2)
        z_npts = (round(z_dist / spacing)) & (-2)

        return [x_npts, y_npts, z_npts]

    def _calculate_rmsd(self, rmsd_path, pose_dir_path, mc=False):
        if not mc:
            self.logger.info('\n#==============================================================================#')
            self.logger.info("CALCULATING RMSD VALUES FOR EACH POSE WITH RESPECT TO THE STARTING XYZ FILE\n")

        print_list = []
        rmsd_list = []
        rmsd_func = Path(os.environ['ROOT_DIR']) / 'metal_dock' / 'calculate_rmsd.py'

        for pose in range(1, self.par.num_poses + 1):
            if not mc:
                pose_path = pose_dir_path / f'pose_{pose}'
            else:
                pose_path = pose_dir_path

            pose_path = pose_path / f'{self.par.name_ligand}_{pose}.xyz'
            command = os.environ['PYTHON_3']+f' {rmsd_func} {rmsd_path} {pose_path} -nh --reorder --rotation none --translation none'
            rmsd_non_rotate = float(subprocess.getoutput([command]))
            rmsd = rmsd_non_rotate

            rmsd_list.append(rmsd)
            print_list.append(f"RMSD for Conformation {pose:>3} = {rmsd:>8.4f}")


        avg_output = np.mean(rmsd_list)
        stdv_rmsd = np.std(rmsd_list)
        var_rmsd = np.var(rmsd_list)
        print_list.append(f'Average RMSD              = {avg_output:8.4f}')
        print_list.append(f'Standard Deviation RMSD   = {stdv_rmsd:8.4f}')
        print_list.append(f"Variance RMSD             = {var_rmsd:8.4f}\n")
        return print_list, rmsd_list, avg_output, stdv_rmsd, var_rmsd


class DockingMC(Docking):
    def __init__(self, par, metal_complex, protein, xyz_path):
        self.xyz_path = xyz_path
        super().__init__(par, metal_complex=None, protein=None)

    def run_mc(self, dock_dir_path):
        """
        Run a MetalDock run for the Monte Carlo optimization protocol.
        """
        gpf_path = dock_dir_path / f'{self.par.name_ligand}.gpf'
        dpf_path = dock_dir_path / f'{self.par.name_ligand}_{self.par.name_protein}.dpf'

        ligand_pdbqt_path = dock_dir_path / f'{self.par.name_ligand}.pdbqt'
        receptor_pdbqt_path = dock_dir_path / f'{self.par.name_protein}.pdbqt'

        if not gpf_path.exists():
            self._create_gpf_file(ligand_pdbqt_path, receptor_pdbqt_path, gpf_path, self.par.parameter_file)
        else:
            os.remove(gpf_path)
            self._create_gpf_file(ligand_pdbqt_path, receptor_pdbqt_path, gpf_path, self.par.parameter_file)

        if not dpf_path.exists():
            self._create_dpf_file(dpf_path, gpf_path, self.par.parameter_file)
        else:
            os.remove(dpf_path)
            self._create_dpf_file(dpf_path, gpf_path, self.par.parameter_file)

        dlg_path = self.par.output_dir / self.par.name_ligand / f'{self.par.name_ligand}_{self.par.name_protein}.dlg'
        autogrid_log_file = self.par.output_dir / self.par.name_ligand / 'autogrid.log'
        autodock_log_file = self.par.output_dir / self.par.name_ligand / 'autodock.log'

        os.chdir(dock_dir_path)
        self._run_autogrid(autogrid_log_file, gpf_path)
        self._run_autodock(autodock_log_file, dpf_path)
        self._write_conformations(dlg_path, dock_dir_path)
        self._write_pdbqt_to_xyz(dock_dir_path)

        rmsd_path = dock_dir_path / f'{self.par.name_ligand}_c.xyz'
        # remove all hydrogens from the rmsd path file 
        self._remove_hydrogens(rmsd_path)
        pose_path = dock_dir_path 

        print_list, rmsd_list, avg_rmsd, stdv_rmsd, var_rmsd = self._calculate_rmsd(rmsd_path, pose_path, mc=True)
        return print_list, rmsd_list, avg_rmsd, stdv_rmsd, var_rmsd

    def _get_box_centre(self):
        """
        Get the box centre for Monte Carlo docking.

        Args:
            metal_symbol (str): The symbol of the metal to get the centre for.

        Returns:
            list: The coordinates of the box centre.
        """
        with open(self.xyz_path, 'r') as f:
            for line in f:
                if self.par.metal_symbol in line:
                    dock_x, dock_y, dock_z = map(float, line.split()[1:4])
                    break
        return [dock_x, dock_y, dock_z]

    def _write_pdbqt_to_xyz(self, dock_dir_path):
        """
        This function writes the pdbqt file to an xyz file without the hydrogens
        """
        for pose in range(1, self.par.num_poses + 1):
            pdqt_path = dock_dir_path / f'{self.par.name_ligand}_{pose}.pdbqt'
            xyz_path = dock_dir_path / f'{self.par.name_ligand}_{pose}.xyz'
            # just write the pdbqt file to an xyz file without the hydrogens
            pdqt_file = open(pdqt_path, 'r')
            xyz_file = open(xyz_path, 'w')
            
            xyz_lines = []
            n_atoms = 0
            for line in pdqt_file:
                if "ATOM" in line and "DD" not in line:
                    xyz_lines.append(f'{line.strip().split()[2]:>2}  {float(line.strip().split()[6]):>8.3f}  {float(line.strip().split()[7]):>8.3f} {float(line.strip().split()[8]):>8.3f}\n')
                    n_atoms += 1
                else:
                    continue

            # insert at the beginning of the xyz_lines the number of atoms
            xyz_lines.insert(0, f'{n_atoms}\n\n')

            # write the xyz_lines to the xyz file
            xyz_file.writelines(xyz_lines)

            pdqt_file.close()
            pdqt_file.close()
            xyz_file.close()

    def _remove_hydrogens(self, xyz_file):
        """
        Remove all hydrogens from the xyz file
        """
        with open(xyz_file, 'r') as f:
            lines = f.readlines()
            lines = [line for line in lines if line.strip().split()[0] != 'H']
        with open(xyz_file, 'w') as f:
            f.writelines(lines)
