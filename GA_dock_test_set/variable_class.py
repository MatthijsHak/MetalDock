class lig_par_dock:

    def __init__(self, metal_symbol, metal_cap, parameter_file, box_size, scale_factor, docking_genetic_algorithm, docking_simulated_annealing,
                    random_position, time_step, max_initial_energy, initial_annealing_temperature, temp_reduction_factor, number_of_runs, max_cycles,
                         r_OA, e_OA, r_SA, e_SA, r_HD, e_HD, r_NA, e_NA, r_N, e_N, r_M, e_M):

        
        self.metal_symbol = metal_symbol
        self.metal_cap = metal_cap
        
        self.parameter_file = parameter_file
        self.box_size = box_size
        self.scale_factor = scale_factor

        self.docking_simulated_annealing = docking_simulated_annealing
        self.docking_genetic_algorithm = docking_genetic_algorithm

        self.random_position = random_position

        if self.docking_simulated_annealing != None:
            self.time_step = time_step
            self.max_initial_energy  = max_initial_energy
            self.initial_annealing_temperature    = initial_annealing_temperature
            self.temp_reduction_factor  = temp_reduction_factor
            self.number_of_runs  = number_of_runs
            self.max_cycles = max_cycles

        if self.docking_genetic_algorithm != None:
            self.docking_genetic_algorithm = docking_genetic_algorithm


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
