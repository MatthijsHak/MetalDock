import sys, os
import re
import configparser

from docking import standard_set 


config = configparser.ConfigParser(interpolation=None)
config['DEFAULT']       =   { "method"                  :       'dock',
                              "metal_symbol"            :           '',
                              "parameter_file"          :'metal_dock.dat',
                              "ncpu"                    :           '1'}

config['PROTEIN']       =   { "pdb_file"                :           '',
                              "pH"                      :        '7.4',
                              "clean_pdb"               :        'True'}

config['METAL_COMPLEX'] =   { "geom_opt"                :       'False',
                              "xyz_file"                :        'None',
                              "charge"                  :           '0',
                              "spin"                    :           '0',
                              "vacant_site"             :         'True',}

config['QM']            =   { "engine"                  :            '',

                              # ADF & Gaussian input keywords
                              "basis_set"               :            '',
                              "functional_type"         :            '',
                              "functional"              :            '',
                              "dispersion"              :            '',
                              "solvent"                 :            '',

                              # ORCA input keywords
                              "orcasimpleinput"         :            '',
                              "orcablocks"              :            ''}

config['DOCKING']       =   { "standard"                :        'True',
                              "rmsd"                    :       'False',
                              "dock_x"                  :         '0.0',
                              "dock_y"                  :        ' 0.0',
                              "dock_z"                  :        ' 0.0',
                              "box_size"                :            '',
                              "scale_factor"            :            '',
                              "random_pos"              :        'True',
                              "num_poses"               :          '10',
                              "e_NA"                    :         '5.0',
                              "e_OA"                    :         '5.0',
                              "e_SA"                    :         '5.0',
                              "e_HD"                    :         '5.0',
                              "ga_dock"                 :       'False',
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

config['MC']            =   { "mc_steps"                 :       '250'}


class Parser:
  
  def __init__(self, input_file):
    config.read(input_file)

    # Convert to correct datatype if necessary
    # [DEFAULT] #
    self.method                   = config['DEFAULT']['method']
    self.metal_symbol             = config['DEFAULT']['metal_symbol']

    self.parameter_file           = config['DEFAULT']['parameter_file']
    self.ncpu                     = int(config['DEFAULT']['ncpu'])

    self.atom_types_included()

    # [PROTEIN] # 
    self.pdb_file                 = config['PROTEIN']['pdb_file']
    self.name_protein             = config['PROTEIN']['pdb_file'][:-4]
    self.pH                       = float(config['PROTEIN']['pH'])
    self.clean_pdb                = config['PROTEIN'].getboolean('clean_pdb')

    # [METAL_COMPLEX] #
    self.geom_opt                 = config['METAL_COMPLEX'].getboolean('geom_opt')
    self.xyz_file                 = config['METAL_COMPLEX']['xyz_file']
    self.name_ligand              = config['METAL_COMPLEX']['xyz_file'][:-4]
    self.charge                   = int(config['METAL_COMPLEX']['charge'])
    self.spin                     = int(config['METAL_COMPLEX']['spin'])
    self.vacant_site              = config['METAL_COMPLEX'].getboolean('vacant_site')

    # [QM] # 
    self.engine                   = str(config['QM']['engine'])
    self.basis_set                = config['QM']['basis_set']
    self.functional_type          = config['QM']['functional_type']
    self.functional               = config['QM']['functional']
    self.dispersion               = config['QM']['dispersion']
    self.solvent                  = config['QM']['solvent']
    self.orcasimpleinput          = config['QM']['orcasimpleinput']
    self.orcablocks               = config['QM']['orcablocks']

    if self.engine == 'ORCA' and self.orcasimpleinput == '':
      print('ORCA engine selected but no ORCA input found')
      print('Please insert QM specifications with the orcasimpleinput & orcablocks keywords in the .ini file')
      sys.exit()

    # [DOCKING] #
    self.standard                 = config['DOCKING'].getboolean('standard')
    if self.standard == False:
      self.parameter_set            = [self.e_NA, self.e_OA, self.e_SA, self.e_HD]

    self.rmsd                     = config['DOCKING'].getboolean('rmsd')
    self.dock_x                   = float(config['DOCKING']['dock_x'])
    self.dock_y                   = float(config['DOCKING']['dock_y'])
    self.dock_z                   = float(config['DOCKING']['dock_z'])
    
    try: 
      self.box_size               = int(config['DOCKING']['box_size'])
    except ValueError:
      self.box_size               = 0 
    
    try:
      self.scale_factor           = float(config['DOCKING']['scale_factor'])
    except ValueError:
      self.scale_factor           = 0 

    self.random_pos               = config['DOCKING'].getboolean('random_pos')
    self.num_poses                = int(config['DOCKING']['num_poses'])
    self.e_NA                     = float(config['DOCKING']['e_NA'])
    self.e_OA                     = float(config['DOCKING']['e_OA'])
    self.e_SA                     = float(config['DOCKING']['e_SA'])
    self.e_HD                     = float(config['DOCKING']['e_HD'])

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

    # [MC] #
    self.mc_steps                 = int(config['MC']['mc_steps'])

    if self.ga_dock == False and self.sa_dock == False:
        print('At least ga_dock or sa_dock must be set to True for MetalDock to run properly')
        sys.exit(1)

  def atom_types_included(self):
    if self.parameter_file == 'metal_dock.dat':
      param_file = os.environ['ROOT_DIR']+'/metal_dock.dat'
    else:
      param_file = self.parameter_file

    # Open the text file for reading
    with open(param_file, 'r') as file:
        # Initialize an empty list to store the symbols
        symbols_list = []

        # Iterate through each line in the file
        for line in file:
            # Use regular expressions to find matches and extract symbols
            matches = re.findall(r'atom_par\s+([A-Za-z]+)', line)
            
            # Extend the symbols list with the matches (if any)
            symbols_list.extend(matches)

    # Remove following atom symbols:
    rmv_list = ['H','HD','HS','C','A','N','NA','NS','OA','OS','F','MG','P','SA','S','Cl','CL','CA','MN','FE','ZN','BR','I','Z','G','GA','J','Q','DD']  

    for rmv in rmv_list:
      if rmv in symbols_list:
        symbols_list.remove(rmv)

    # Print the list of symbols
    if self.metal_symbol not in symbols_list and self.parameter_file == 'metal_dock.dat':
      print('The metal symbol you have chosen is currently not supported by MetalDock')
      print('The following atom types are present in the parameter file: ')
      print(' '.join(symbols_list))
      sys.exit()
    elif self.metal_symbol not in symbols_list and self.parameter_file != 'metal_dock.dat':
      print('The metal symbol you have chosen is currently not supported by your selected parameter file')
      print('The following atom types are present in the parameter file: ')
      print(' '.join(symbols_list))
      print('Please choose a different parameter file or add the missing atom types to the parameter file')
      sys.exit()