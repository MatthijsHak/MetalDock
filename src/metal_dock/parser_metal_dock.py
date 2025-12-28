import sys
import os
import re
import configparser
from pathlib import Path
from src.metal_dock.logger import MetalDockLogger
from src.metal_dock.environment_variables import ROOT_DIR

# Define standard default values
STANDARD_DEFAULTS = {
    'method': '',
    'metal_symbol': 'Ru',
    'parameter_file': f'{os.environ["ROOT_DIR"]}/metal_dock/metal_dock.dat',
    'ncpu': 1,
    'memory': 1500,
    'pdb_file': 'protein.pdb',
    'pH': 7.4,
    'clean_pdb': True,
    'geom_opt': True,
    'xyz_file': 'metal_complex.xyz',
    'charge': 0,
    'spin': 0,
    'vacant_site': True,
    'engine': 'ADF',
    'basis_set': 'TZP',
    'functional_type': 'GGA',
    'functional': 'PBE',
    'dispersion': '',
    'relativity': '',
    'solvent': '',
    'orcasimpleinput': 'PBE def2-TZVP CPCM(Water)',
    'orcablocks': '',
    'ini_parameters': False,
    'rmsd': False,
    'dock_x': None,
    'dock_y': None,
    'dock_z': None,
    'box_size': 20,
    'scale_factor': 0.0,
    'random_pos': True,
    'num_poses': 10,
    'e_NA': 5.0,
    'e_OA': 5.0,
    'e_SA': 5.0,
    'e_HD': 5.0,
    'ga_dock': True,
    'ga_dock_pop_size': 150,
    'ga_dock_num_evals': 2500000,
    'ga_dock_num_generations': 27000,
    'ga_dock_elitism': 1,
    'ga_dock_mutation_rate': 0.02,
    'ga_dock_crossover_rate': 0.80,
    'ga_dock_window_size': 10,
    'sa_dock': False,
    'temp_reduction_factor': 0.90,
    'number_of_runs': 50,
    'max_cycles': 50,
    'mc_steps': 250,
    'score_only': False
}

             #         e_NA, e_OA, e_SA,  e_HD,
standard_set = {'V' : [ 4.696,	6.825,	5.658,	3.984],
                'CR': [ 6.371,	1.998,	0.144,	3.625],
                'CO': [ 5.280,	0.050,	6.673,	5.929],
                'NI': [ 0.630,	2.732,	4.462,	2.820],
                'CU': [ 4.696,	1.277,	6.791,	1.114],
                'MO': [ 1.330,	0.014,	0.168,	5.620],
                'RU': [ 6.936,	2.796,	4.295,	6.357],
                'RH': [ 5.559,	2.056,	0.573,	5.471],
                'PD': [ 4.688,	0.845,	5.574,	3.159],
                'RE': [ 6.738,	0.645,	3.309,	4.502],
                'OS': [ 5.958,	0.135,	4.102,	6.589],
                'PT': [ 6.532,	2.020,	6.332,	1.844],
            }


class ParserBase:
    def __init__(self, config_file):
        self.config = configparser.ConfigParser(interpolation=None)
        self.config.read(config_file)
        self.logger = MetalDockLogger()

    def get_config_value(self, section, key, default=None, cast_type=None):
        """
        Helper function to get the value from the config file.
        """
        value = self.config[section].get(key, default)
        if cast_type and value is not None:
            if cast_type is bool:
                return str(value).strip().lower() in ("true", "1", "yes")
            return cast_type(value)
        return value

    def find_metal_symbol(self, xyz_file, metal_symbol):
        with open(xyz_file, 'r') as file:
            for line in file.readlines()[2:]:
                if line.strip() == '':
                    break
                if line.strip().split()[0] == metal_symbol:
                    return True
        return False

    def calculate_heavy_atoms(self, xyz_file):
        """
        Calculate the number of heavy atoms in the xyz file.
        """
        with open(xyz_file, 'r') as file:
            return sum(1 for line in file.readlines()[2:] if line.strip() and line.strip().split()[0] != 'H')

    def atom_types_included(self):
        """
        Verify that the atom types are present in the parameter file.
        """
        param_file = self.parameter_file if self.parameter_file != 'metal_dock.dat' else os.path.join(os.environ['ROOT_DIR'], 'metal_dock/metal_dock.dat')
        with open(param_file, 'r') as file:
            symbols_list = [match for line in file for match in re.findall(r'atom_par\s+([A-Za-z]+)', line)]

        rmv_list = ['H', 'HD', 'HS', 'C', 'A', 'N', 'NA', 'NS', 'OA', 'OS', 'F', 'MG', 'P', 'SA', 'S', 'Cl', 'CL', 'CA', 'MN', 'FE', 'ZN', 'BR', 'I', 'Z', 'G', 'GA', 'J', 'Q', 'DD']
        symbols_list = [symbol for symbol in symbols_list if symbol not in rmv_list]

        if self.metal_symbol not in symbols_list:
            self.logger.info(f'THE METAL SYMBOL YOU HAVE CHOSEN IS NOT SUPPORTED BY THE SELECTED PARAMETER FILE')
            self.logger.info('THE FOLLOWING ATOM TYPES ARE PRESENT IN THE PARAMETER FILE: ')
            self.logger.info(' '.join(symbols_list))
            if self.parameter_file != 'metal_dock.dat':
                self.logger.info('PLEASE CHOOSE A DIFFERENT PARAMETER FILE OR ADD THE MISSING ATOM TYPES TO THE PARAMETER FILE')
            sys.exit()

    def _parse_box_size(self):
        """
        Parse the box size for the DockParser class.
        """
        box_size = self.config['DOCKING'].get('box_size', str(STANDARD_DEFAULTS['box_size'])).strip().split(',')
        if len(box_size) == 1:
            return [float(box_size[0])] * 3
        elif len(box_size) == 2:
            self.logger.info('ONLY TWO BOX SIZES GIVEN, SETTING THIRD BOX SIZE TO 20\n')
            return [float(size) for size in box_size] + [20]
        return [float(size) for size in box_size]

    def _validate_docking_algorithm(self):
        """
        Validate the docking algorithm.
        """
        if self.ga_dock and self.sa_dock:
            self.logger.info('Only one docking algorithm can be chosen, set either ga_dock or sa_dock to False')
            sys.exit()

        if self.ga_dock:
            self.dock_algorithm = [self.ga_dock_pop_size, self.ga_dock_num_evals, self.ga_dock_num_generations,
                                   self.ga_dock_elitism, self.ga_dock_mutation_rate, self.ga_dock_crossover_rate,
                                   self.ga_dock_window_size]

        if self.sa_dock:
            self.dock_algorithm = [self.temp_reduction_factor, self.number_of_runs, self.max_cycles]

        if self.score_only and (self.ga_dock or self.sa_dock):
            self.logger.info('ga_dock and/or sa_dock must be set to False when performing score_only')
            sys.exit()

        if not self.ga_dock and not self.sa_dock and not self.score_only:
            self.logger.info('At least ga_dock or sa_dock or score_only must be set to True for MetalDock to run properly')
            sys.exit()


class MCParser(ParserBase):
    def __init__(self, config_file):
        super().__init__(config_file)
        self._initialize_parameters()

    def _initialize_parameters(self):
        """
        Initialize the parameters for the MCParser class.
        """
        # [DEFAULT] #
        self.metal_symbol = self.get_config_value('DEFAULT', 'metal_symbol', STANDARD_DEFAULTS['metal_symbol']).strip()
        self.ncpu = self.get_config_value('DEFAULT', 'ncpu', STANDARD_DEFAULTS['ncpu'], int)
        self.memory = self.get_config_value('DEFAULT', 'memory', STANDARD_DEFAULTS['memory'], int)

        # [MC] #
        self.mc_steps = self.get_config_value('MC', 'mc_steps', STANDARD_DEFAULTS['mc_steps'], int)
        self.parameter_file = Path(self.get_config_value('MC', 'parameter_file', STANDARD_DEFAULTS['parameter_file']).strip()).resolve()

        # [DOCKING] #
        self.dock_x = self.get_config_value('DOCKING', 'dock_x', STANDARD_DEFAULTS['dock_x'], float)
        self.dock_y = self.get_config_value('DOCKING', 'dock_y', STANDARD_DEFAULTS['dock_y'], float)
        self.dock_z = self.get_config_value('DOCKING', 'dock_z', STANDARD_DEFAULTS['dock_z'], float)

        self.box_size = self._parse_box_size()
        self.scale_factor = self.get_config_value('DOCKING', 'scale_factor', STANDARD_DEFAULTS['scale_factor'], float)
        self.random_pos = self.get_config_value('DOCKING', 'random_pos', STANDARD_DEFAULTS['random_pos'], bool)
        self.num_poses = self.get_config_value('DOCKING', 'num_poses', STANDARD_DEFAULTS['num_poses'], int)
        self.e_NA = self.get_config_value('DOCKING', 'e_NA', STANDARD_DEFAULTS['e_NA'], float)
        self.e_OA = self.get_config_value('DOCKING', 'e_OA', STANDARD_DEFAULTS['e_OA'], float)
        self.e_SA = self.get_config_value('DOCKING', 'e_SA', STANDARD_DEFAULTS['e_SA'], float)
        self.e_HD = self.get_config_value('DOCKING', 'e_HD', STANDARD_DEFAULTS['e_HD'], float)

        self.ga_dock = self.get_config_value('DOCKING', 'ga_dock', STANDARD_DEFAULTS['ga_dock'], bool)
        self.ga_dock_pop_size = self.get_config_value('DOCKING', 'ga_dock_pop_size', STANDARD_DEFAULTS['ga_dock_pop_size'], int)
        self.ga_dock_num_evals = self.get_config_value('DOCKING', 'ga_dock_num_evals', STANDARD_DEFAULTS['ga_dock_num_evals'], int)
        self.ga_dock_num_generations = self.get_config_value('DOCKING', 'ga_dock_num_generations', STANDARD_DEFAULTS['ga_dock_num_generations'], int)
        self.ga_dock_elitism = self.get_config_value('DOCKING', 'ga_dock_elitism', STANDARD_DEFAULTS['ga_dock_elitism'], int)
        self.ga_dock_mutation_rate = self.get_config_value('DOCKING', 'ga_dock_mutation_rate', STANDARD_DEFAULTS['ga_dock_mutation_rate'], float)
        self.ga_dock_crossover_rate = self.get_config_value('DOCKING', 'ga_dock_crossover_rate', STANDARD_DEFAULTS['ga_dock_crossover_rate'], float)
        self.ga_dock_window_size = self.get_config_value('DOCKING', 'ga_dock_window_size', STANDARD_DEFAULTS['ga_dock_window_size'], int)

        self.sa_dock = self.get_config_value('DOCKING', 'sa_dock', STANDARD_DEFAULTS['sa_dock'], bool)
        self.temp_reduction_factor = self.get_config_value('DOCKING', 'temp_reduction_factor', STANDARD_DEFAULTS['temp_reduction_factor'], float)
        self.number_of_runs = self.get_config_value('DOCKING', 'number_of_runs', STANDARD_DEFAULTS['number_of_runs'], int)
        self.max_cycles = self.get_config_value('DOCKING', 'max_cycles', STANDARD_DEFAULTS['max_cycles'], int)

        self._validate_docking_algorithm()

class DockParser(ParserBase):
    def __init__(self, config_file):
        super().__init__(config_file)
        self._initialize_parameters()

    def _initialize_parameters(self):
        """
        Initialize the parameters for the DockParser class.
        """
        # [DEFAULT] #
        self.metal_symbol = self.get_config_value('DEFAULT', 'metal_symbol', STANDARD_DEFAULTS['metal_symbol']).strip()
        self.parameter_file = Path(self.get_config_value('DEFAULT', 'parameter_file', STANDARD_DEFAULTS['parameter_file']).strip()).resolve()
        self.ncpu = self.get_config_value('DEFAULT', 'ncpu', STANDARD_DEFAULTS['ncpu'], int)
        self.memory = self.get_config_value('DEFAULT', 'memory', STANDARD_DEFAULTS['memory'], int)

        # [PROTEIN] #
        self.pdb_file = self.get_config_value('PROTEIN', 'pdb_file', STANDARD_DEFAULTS['pdb_file']).strip()
        self.name_protein = self.pdb_file.split('/')[-1][:-4]
        self.pH = self.get_config_value('PROTEIN', 'pH', STANDARD_DEFAULTS['pH'], float)
        self.clean_pdb = self.get_config_value('PROTEIN', 'clean_pdb', STANDARD_DEFAULTS['clean_pdb'], bool)

        # [METAL_COMPLEX] #
        self.geom_opt = self.get_config_value('METAL_COMPLEX', 'geom_opt', STANDARD_DEFAULTS['geom_opt'], bool)
        self.xyz_file = self.get_config_value('METAL_COMPLEX', 'xyz_file', STANDARD_DEFAULTS['xyz_file']).strip()
        self.name_ligand = self.xyz_file.split('/')[-1][:-4]
        self.charge = self.get_config_value('METAL_COMPLEX', 'charge', STANDARD_DEFAULTS['charge'], int)
        self.spin = self.get_config_value('METAL_COMPLEX', 'spin', STANDARD_DEFAULTS['spin'], float)
        self.vacant_site = self.get_config_value('METAL_COMPLEX', 'vacant_site', STANDARD_DEFAULTS['vacant_site'], bool)

        # [QM] #
        self.engine = self.get_config_value('QM', 'engine', STANDARD_DEFAULTS['engine']).strip()
        self.basis_set = self.get_config_value('QM', 'basis_set', STANDARD_DEFAULTS['basis_set']).strip()
        self.functional_type = self.get_config_value('QM', 'functional_type', STANDARD_DEFAULTS['functional_type']).strip()
        self.functional = self.get_config_value('QM', 'functional', STANDARD_DEFAULTS['functional']).strip()
        self.dispersion = self.get_config_value('QM', 'dispersion', STANDARD_DEFAULTS['dispersion']).strip()
        self.relativity = self.get_config_value('QM', 'relativity', STANDARD_DEFAULTS['relativity']).strip()
        self.solvent = self.get_config_value('QM', 'solvent', STANDARD_DEFAULTS['solvent']).strip()
        self.orcasimpleinput = self.get_config_value('QM', 'orcasimpleinput', STANDARD_DEFAULTS['orcasimpleinput']).strip()
        self.orcablocks = self.get_config_value('QM', 'orcablocks', STANDARD_DEFAULTS['orcablocks']).strip()

        if self.engine == 'ORCA' and not self.orcasimpleinput:
            self.logger.info('ORCA engine selected but no ORCA input found')
            self.logger.info('Please insert QM specifications with the orcasimpleinput & orcablocks keywords in the .ini file')
            sys.exit()

        # [DOCKING] #
        self.ini_parameters = self.get_config_value('DOCKING', 'ini_parameters', STANDARD_DEFAULTS['ini_parameters'], bool)
        self.internal_param = self.metal_symbol.upper() in ['FE', 'ZN', 'MN']
        self.parameter_set = self._get_parameter_set()

        self.rmsd = self.get_config_value('DOCKING', 'rmsd', STANDARD_DEFAULTS['rmsd'], bool)
        self.dock_x = self.get_config_value('DOCKING', 'dock_x', STANDARD_DEFAULTS['dock_x'], float)
        self.dock_y = self.get_config_value('DOCKING', 'dock_y', STANDARD_DEFAULTS['dock_y'], float)
        self.dock_z = self.get_config_value('DOCKING', 'dock_z', STANDARD_DEFAULTS['dock_z'], float)

        self.box_size = self._parse_box_size()
        self.scale_factor = self.get_config_value('DOCKING', 'scale_factor', STANDARD_DEFAULTS['scale_factor'], float)
        self.random_pos = self.get_config_value('DOCKING', 'random_pos', STANDARD_DEFAULTS['random_pos'], bool)
        self.num_poses = self.get_config_value('DOCKING', 'num_poses', STANDARD_DEFAULTS['num_poses'], int)
        self.e_NA = self.get_config_value('DOCKING', 'e_NA', STANDARD_DEFAULTS['e_NA'], float)
        self.e_OA = self.get_config_value('DOCKING', 'e_OA', STANDARD_DEFAULTS['e_OA'], float)
        self.e_SA = self.get_config_value('DOCKING', 'e_SA', STANDARD_DEFAULTS['e_SA'], float)
        self.e_HD = self.get_config_value('DOCKING', 'e_HD', STANDARD_DEFAULTS['e_HD'], float)

        self.ga_dock = self.get_config_value('DOCKING', 'ga_dock', STANDARD_DEFAULTS['ga_dock'], bool)
        self.ga_dock_pop_size = self.get_config_value('DOCKING', 'ga_dock_pop_size', STANDARD_DEFAULTS['ga_dock_pop_size'], int)
        self.ga_dock_num_evals = self.get_config_value('DOCKING', 'ga_dock_num_evals', STANDARD_DEFAULTS['ga_dock_num_evals'], int)
        self.ga_dock_num_generations = self.get_config_value('DOCKING', 'ga_dock_num_generations', STANDARD_DEFAULTS['ga_dock_num_generations'], int)
        self.ga_dock_elitism = self.get_config_value('DOCKING', 'ga_dock_elitism', STANDARD_DEFAULTS['ga_dock_elitism'], int)
        self.ga_dock_mutation_rate = self.get_config_value('DOCKING', 'ga_dock_mutation_rate', STANDARD_DEFAULTS['ga_dock_mutation_rate'], float)
        self.ga_dock_crossover_rate = self.get_config_value('DOCKING', 'ga_dock_crossover_rate', STANDARD_DEFAULTS['ga_dock_crossover_rate'], float)
        self.ga_dock_window_size = self.get_config_value('DOCKING', 'ga_dock_window_size', STANDARD_DEFAULTS['ga_dock_window_size'], int)

        self.sa_dock = self.get_config_value('DOCKING', 'sa_dock', STANDARD_DEFAULTS['sa_dock'], bool)
        self.temp_reduction_factor = self.get_config_value('DOCKING', 'temp_reduction_factor', STANDARD_DEFAULTS['temp_reduction_factor'], float)
        self.number_of_runs = self.get_config_value('DOCKING', 'number_of_runs', STANDARD_DEFAULTS['number_of_runs'], int)
        self.max_cycles = self.get_config_value('DOCKING', 'max_cycles', STANDARD_DEFAULTS['max_cycles'], int)

        self.score_only = self.get_config_value('DOCKING', 'score_only', STANDARD_DEFAULTS['score_only'], bool)

        self._validate_docking_algorithm()

        # heavy atoms
        self.n_heavy_atoms = self.calculate_heavy_atoms(self.xyz_file)

        # check if metal symbol is present in the xyz file
        if not self.find_metal_symbol(self.xyz_file, self.metal_symbol):
            self.logger.info('The metal symbol you have chosen is not present in the xyz file')
            self.logger.info('Please choose a different metal symbol or add the missing metal symbol to the xyz file')
            sys.exit()

    def _get_parameter_set(self):
        """
        Get the parameter set for the DockParser class.
        """
        if self.internal_param:
            return None
        if self.ini_parameters:
            return [
                self.get_config_value('DOCKING', 'e_NA', STANDARD_DEFAULTS['e_NA'], float),
                self.get_config_value('DOCKING', 'e_OA', STANDARD_DEFAULTS['e_OA'], float),
                self.get_config_value('DOCKING', 'e_SA', STANDARD_DEFAULTS['e_SA'], float),
                self.get_config_value('DOCKING', 'e_HD', STANDARD_DEFAULTS['e_HD'], float)
            ]
        return standard_set.get(self.metal_symbol.upper(), [5.0, 5.0, 5.0, 5.0])
