class lig_par_dock_SA:

    def __init__(self, standard, rmsd, dock_x, dock_y, dock_z, pdb_file_protein, name_protein, pH, xyz_file_ligand, name_ligand, parameter_file, scale_factor, box_size, metal_symbol, metal_cap, charge_ligand, 
                        spin_ligand,  basis_set, gga_functional, hybrid_functional, dispersion_correction, docking_simulated_annealing, docking_genetic_algorithm, random_position, time_step, max_initial_energy, initial_annealing_temperature, temp_reduction_factor, 
                            number_of_runs, max_cycles, r_OA, e_OA, r_SA, e_SA, r_HD, e_HD, r_NA, e_NA, r_N, e_N, r_M, e_M):

        self.standard = standard
        self.rmsd = rmsd

        self.dock_x = dock_x
        self.dock_y = dock_y
        self.dock_z = dock_z

        self.pdb_file_protein = pdb_file_protein
        self.name_protein = name_protein
        self.pH = pH

        self.xyz_file_ligand = xyz_file_ligand
        self.name_ligand = name_ligand

        self.parameter_file = parameter_file
        self.scale_factor = scale_factor
        self.box_size = box_size

        self.metal_symbol = metal_symbol
        self.metal_cap = metal_cap
        self.charge_ligand = charge_ligand
        self.spin_ligand = spin_ligand

        self.basis_set = basis_set
        self.gga_functional = gga_functional
        self.hybrid_functional = hybrid_functional
        self.dispersion_correction = dispersion_correction

        self.random_position = random_position
        self.docking_simulated_annealing = docking_simulated_annealing
        self.docking_genetic_algorithm = docking_genetic_algorithm

        self.time_step = time_step
        self.max_initial_energy  = max_initial_energy
        self.initial_annealing_temperature    = initial_annealing_temperature
        self.temp_reduction_factor  = temp_reduction_factor
        self.number_of_runs  = number_of_runs
        self.max_cycles = max_cycles

        self.r_OA = r_OA
        self.e_OA = e_OA

        self.r_SA = r_SA
        self.e_SA = e_SA

        self.r_HD = r_HD
        self.e_HD = e_HD

        self.r_NA = r_NA
        self.e_NA = e_NA

        self.r_N = r_N
        self.e_N = e_N

        self.r_M = r_M
        self.e_M = e_M


class lig_par_dock_GA:

    def __init__(self, standard, rmsd, dock_x, dock_y, dock_z, pdb_file_protein, name_protein, pH, xyz_file_ligand, name_ligand, parameter_file, scale_factor, box_size, metal_symbol, metal_cap, charge_ligand, 
                        spin_ligand,  basis_set, gga_functional, hybrid_functional, dispersion_correction, docking_simulated_annealing, docking_genetic_algorithm, random_position, ga_pop_size, ga_num_evals, ga_num_generation, ga_elitism, 
                            ga_mutation_rate, ga_crossover_rate, ga_window_size, r_OA, e_OA, r_SA, e_SA, r_HD, e_HD, r_NA, e_NA, r_N, e_N, r_M, e_M):

        self.standard = standard
        self.rmsd = rmsd

        self.dock_x = dock_x
        self.dock_y = dock_y
        self.dock_z = dock_z

        self.pdb_file_protein = pdb_file_protein
        self.name_protein = name_protein
        self.pH = pH

        self.xyz_file_ligand = xyz_file_ligand
        self.name_ligand = name_ligand

        self.parameter_file = parameter_file
        self.scale_factor = scale_factor
        self.box_size = box_size

        self.metal_symbol = metal_symbol
        self.metal_cap = metal_cap
        self.charge_ligand = charge_ligand
        self.spin_ligand = spin_ligand

        self.basis_set = basis_set
        self.gga_functional = gga_functional
        self.hybrid_functional = hybrid_functional
        self.dispersion_correction = dispersion_correction

        self.random_position = random_position
        self.docking_simulated_annealing = docking_simulated_annealing
        self.docking_genetic_algorithm = docking_genetic_algorithm
        
        self.ga_pop_size = ga_pop_size
        self.ga_num_evals = ga_num_evals
        self.ga_num_generation  = ga_num_generation
        self.ga_elitism = ga_elitism
        self.ga_mutation_rate = ga_mutation_rate
        self.ga_crossover_rate = ga_crossover_rate
        self.ga_window_size = ga_window_size


        self.r_OA = r_OA
        self.e_OA = e_OA

        self.r_SA = r_SA
        self.e_SA = e_SA

        self.r_HD = r_HD
        self.e_HD = e_HD

        self.r_NA = r_NA
        self.e_NA = e_NA

        self.r_N = r_N
        self.e_N = e_N

        self.r_M = r_M
        self.e_M = e_M
