import sys
import configparser


"""
All parameters that can be specified for the following three options:

(1) Docking with metal ligand
(2) Training the GA for a speficic metal 
(3) Testing the GA for a specific metal

# Atom symbol of metal. 
metal_symbol = string

# Parameter file used for docking / GA run. 
parameter_file = string 

# If True docking procedure with standard parameters.
std = bool

# If True calculate the rmsd of the docking procedure with the initial xyz file as reference.
rmsd = bool

#### Protein settings ####
# The path to the pdb file of the protein.
pdb_file = string

# The pH at which the protein structure was determined
pH = float

#### Ligand settings ####
# The path to the xyz file of the ligand
xyz = string

# The charge of the ligand
charge = int

# The spin of the ligand
spin = int

#### QM settings ####
# The quantum software package used. Can be adf, gaussian or psi4 
engine = string

# Basis set used for QM calculation
basis_set = string

# Functional used for QM calculation
functional = string

# Dispersion corrections for QM calculation
dispersion = string

#### Docking settings ####
# The x-coordinate of the centre of the docking box
dock_x = float 

# The y-coordinate of the centre of the docking box
dock_y = float 

# The z-coordinate of the centre of the docking box
dock_z = float 

# Box size is specifed by the number of x-, y- and z-grid points. Each must be an even integer number and are separated by a comma.
npts = int,int,int 

# Box size is cubic and the length of the axis is determined by the maximum distance between two atoms of the molecule. 
# The box size can then be scaled by multiplying this size with a float
scale_factor = float

# The position of the molecule will be randomized if True.
random_pos = bool

# The number of poses that will be docked.
num_poses = int

#### Docking Parameters ####
# LJ parameters for a metal atom and a oxygen that accepts a hydrogen bond.
r_OA = float
e_OA = float

# LJ parameters for a metal atom and a sulfur that accepts a hdyrogen bond. 
r_SA = float
e_SA = float

# LJ parameters for a metal atom and a hydrogen that donates a hydrogen bond. 
r_HD = float
e_HD = float

# LJ parameters for a metal atom and a nitrogen that accepts a hydrogen bond. 
r_NA = float
e_NA = float

# LJ parameters for a metal atom and a nitrogen that does not accept a hydrogen bond. 
r_N = float
e_N = float

# LJ parameters for a metal atom with itself. 
r_M = float
e_M = float

#### Docking Algorithm ####
  ## Genetic Algorithm ##

# If True a genetic algorithm will be used for the docking procedure.
ga_dock = bool 

# The population size of the genetic algorithm for the docking procedure.
ga_dock_pop_size = int

# The maximum number of energy evluations of the genetic algorithm for the docking procedure.
ga_dock_num_evals = int

# The maximum number of generations of the genetic algorithm for the docking procedure.
ga_dock_num_generations = int 

# The number of parents kept for the next generation of the genetic algorithm for the docking procedure.
ga_dock_elitism = int

# The mutation rate of the genetic algorithm for the docking procedure.
ga_dock_mutation_rate = float

# The crossover rate of the genetic algorithm for the docking procedure.
ga_dock_crossover_rate = float

# The window size of the genetic algorithm for the docking procedure.
ga_dock_window_size = float

  ## Simulated Annealing ##

# If True a simulated annealing algorithm will be used for the docking procedure.
sa_dock = bool 

# The maximum initial energy
max_initial_energy = int

# The annealing reduction factor
temp_reduction_factor = float

# The maximum number of docking runs
number_of_runs = int

# The maximum number of temperature reduction cycles
max_cycles = float

###### Parameters specific only for option (2) - TrainParser ######
# The box size is increased during the traning session if the LJ parameters stop changing.
step_wise = bool

# The maximum number of generations for the genetic algorithm to obtain the LJ parameters 
ga_num_generations = int

# The number of parents mating for the genetic algorithm to obtain the LJ parameters
ga_num_mating = int

# The number of solutions within the population for the genetic algorithm to obtain the LJ parameters
sol_per_pop = int

# The parent selection type
par_type = string

# The number of parentes to keep from one generation to the next
keep_par = int 

# If selection type is k_tournamet then the number of parents that will take part in tournament selection
k_tour = int

# The crossover operation type
crossover_type = string

# The probability that the crossover operation will take place 
cross_prob = float [0.0:1.0]

# The probability that the mutation operation will take place
mut_prob = float [0.0:1.0]

"""
config = configparser.ConfigParser()
config['DEFAULT'] = { "metal_symbol"            :            '',
                      "parameter_file"          :            '' }

config['PROTEIN'] = { "pdb_file"                :           '',
                      "pH"                      :        '7.4',
                      "clean_pdb"               :        'True'}

config['LIGAND'] =  { "geom_opt"                :       'False',
                      "xyz_file"                :        'None',
                      "charge"                  :           '0',
                      "spin"                    :           '0'}

config['QM'] =      { "engine"                  :        'adf',
                      "basis_set"               :         'TZP',
                      "functional_type"         :      'hybrid',
                      "functional"              :       'B3LYP',
                      "dispersion"              :            ''}

config['DOCKING'] = { "std"                     :        'True',
                      "rmsd"                    :       'False',
                      "dock_x"                  :         '0.0',
                      "dock_y"                  :        ' 0.0',
                      "dock_z"                  :        ' 0.0',
                      "box_size"                :          '54',
                      "scale_factor"            :           '0',
                      "random_pos"              :        'True',
                      "num_poses"               :          '10',
                      "r_OA"                    :         '2.0',
                      "e_OA"                    :        '10.0',
                      "r_SA"                    :         '2.0',
                      "e_SA"                    :        '10.0',
                      "r_HD"                    :         '2.0',
                      "e_HD"                    :        '10.0',
                      "r_NA"                    :         '2.0',
                      "e_NA"                    :        '10.0',
                      "r_N"                     :         '2.0',
                      "e_N"                     :        '10.0',
                      "r_M"                     :         '2.0',
                      "e_M"                     :        '10.0',
                      "ga_dock"                 :        'True',
                      "ga_dock_pop_size"        :         '150',
                      "ga_dock_num_evals"       :     '2500000',
                      "ga_dock_num_generations" :       '27000',
                      "ga_dock_elitism"         :           '1',
                      "ga_dock_mutation_rate"   :        '0.02',
                      "ga_dock_crossover_rate"  :        '0.80',
                      "ga_dock_window_size"     :          '10',
                      "sa_dock"                 :       'False',
                      "temp_reduction_factor"   :        '0.90',
                      "number_of_runs"          :          '50',
                      "max_cycles"              :          '50'}

config['GA']  =     { "step_wise"               :        'True',
                      "ga_num_generations"      :          '20',
                      "ga_num_mating"           :           '2',
                      "sol_per_pop"             :          '20',
                      "par_type"                :         'sss',
                      "keep_par"                :           '1',
                      "k_tour"                  :           '0',
                      "crossover_type"          :     'uniform',
                      "cross_prob"              :        '0.50',
                      "mut_prob"                :        '0.30'}


class Parser:
  
  def __init__(self, input_file):
    config.read(input_file)

    # Convert to correct datatype if necessary
    # [DEFAULT] #
    self.metal_symbol             = config['DEFAULT']['metal_symbol']
    self.parameter_file           = config['DEFAULT']['parameter_file']

    # [PROTEIN] # 
    self.pdb_file                 = config['PROTEIN']['pdb_file']
    self.name_protein             = config['PROTEIN']['pdb_file'][:-4]
    self.pH                       = float(config['PROTEIN']['pH'])
    self.clean_pdb                = config['PROTEIN'].getboolean('clean_pdb')

    # [LIGAND] #
    self.geom_opt                 = config['LIGAND'].getboolean('geom_opt')
    self.xyz_file                 = config['LIGAND']['xyz_file']
    self.name_ligand              = config['LIGAND']['xyz_file'][:-4]
    self.charge                   = int(config['LIGAND']['charge'])
    self.spin                     = int(config['LIGAND']['spin'])

    # [QM] # 
    self.engine                   = str(config['QM']['engine'])
    self.basis_set                = config['QM']['basis_set']
    self.functional_type          = config['QM']['functional_type']
    self.functional               = config['QM']['functional']
    self.dispersion               = config['QM']['dispersion']

    # [DOCKING] #
    self.standard                 = config['DOCKING'].getboolean('std')
    self.rmsd                     = config['DOCKING'].getboolean('rmsd')
    self.dock_x                   = float(config['DOCKING']['dock_x'])
    self.dock_y                   = float(config['DOCKING']['dock_y'])
    self.dock_z                   = float(config['DOCKING']['dock_z'])
    self.box_size                 = int(config['DOCKING']['box_size'])
    self.scale_factor             = float(config['DOCKING']['scale_factor'])
    self.random_pos               = config['DOCKING'].getboolean('random_pos')
    self.num_poses                = int(config['DOCKING']['num_poses'])
    self.r_OA                     = float(config['DOCKING']['r_OA'])
    self.e_OA                     = float(config['DOCKING']['e_OA'])
    self.r_SA                     = float(config['DOCKING']['r_SA'])
    self.e_SA                     = float(config['DOCKING']['e_SA'])
    self.r_HD                     = float(config['DOCKING']['r_HD'])
    self.e_HD                     = float(config['DOCKING']['e_HD'])
    self.r_NA                     = float(config['DOCKING']['r_NA'])
    self.e_NA                     = float(config['DOCKING']['e_NA'])
    self.r_N                      = float(config['DOCKING']['r_N'])
    self.e_N                      = float(config['DOCKING']['e_N'])
    self.r_M                      = float(config['DOCKING']['r_M'])
    self.e_M                      = float(config['DOCKING']['e_M'])

    self.parameter_set            = [self.r_OA, self.e_OA, self.r_SA, self.e_SA, self.r_HD, self.e_HD, self.r_NA, self.e_NA, self.r_N, self.e_N, self.r_M, self.e_M]

    self.ga_dock                  = config['DOCKING'].getboolean('ga_dock')
    self.ga_dock_pop_size         = int(config['DOCKING']['ga_dock_pop_size'])
    self.ga_dock_num_evals        = int(config['DOCKING']['ga_dock_num_evals'])   
    self.ga_dock_num_generations  = int(config['DOCKING']['ga_dock_num_generations'])
    self.ga_dock_elitism          = int(config['DOCKING']['ga_dock_elitism'])
    self.ga_dock_mutation_rate    = float(config['DOCKING']['ga_dock_mutation_rate'])
    self.ga_dock_crossover_rate   = float(config['DOCKING']['ga_dock_crossover_rate'])
    self.ga_dock_window_size      = int(config['DOCKING']['ga_dock_window_size'])

    if self.ga_dock == True:
      self.dock_algorithm         = [self.ga_dock_pop_size, self.ga_dock_num_evals, self.ga_dock_num_generations, self.ga_dock_elitism, self.ga_dock_mutation_rate, self.ga_dock_crossover_rate, self.ga_dock_window_size]

    self.sa_dock                  = config['DOCKING'].getboolean('sa_dock')
    self.temp_reduction_factor    = float(config['DOCKING']['temp_reduction_factor'])
    self.number_of_runs           = int(config['DOCKING']['number_of_runs'])
    self.max_cycles               = int(config['DOCKING']['max_cycles'])

    if self.sa_dock == True:
      self.dock_algorithm        = [self.temp_reduction_factor, self.number_of_runs, self.max_cycles]

    # [GA] # 
    self.step_wise               = config['GA'].getboolean('step_wise')
    self.ga_num_generations      = int(config['GA']['ga_num_generations'])
    self.ga_num_mating           = int(config['GA']['ga_num_mating'])
    self.sol_per_pop             = int(config['GA']['sol_per_pop'])
    self.par_type                = config['GA']['par_type']
    self.keep_par                = int(config['GA']['keep_par'])
    self.k_tour                  = int(config['GA']['k_tour'])
    self.crossover_type          = config['GA']['crossover_type']
    self.cross_prob              = float(config['GA']['cross_prob'])
    self.mut_prob                = float(config['GA']['mut_prob'])