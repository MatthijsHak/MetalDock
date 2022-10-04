class variables:

    def __init__(self, parameter_file, box_size,  scale_factor, step_wise, metal_symbol, metal_cap, num_generations, num_parents_mating, sol_per_pop, parent_selection_type, keep_parents, k_tournament, crossover_type, 
                        crossover_prob, mutation_prob, docking_simulated_annealing, docking_genetic_algorithm, random_position, time_step, max_initial_energy, initial_annealing_temperature, temp_reduction_factor, number_of_runs, max_cycles):

        self.parameter_file = parameter_file
        self.box_size = box_size
        self.scale_factor = scale_factor
        self.step_wise = step_wise

        self.metal_symbol = metal_symbol
        self.metal_cap = metal_cap

        self.sol_per_pop = sol_per_pop
        self.num_generations = num_generations
        self.num_parents_mating = num_parents_mating

        self.parent_selection_type = parent_selection_type
        self.keep_parents = keep_parents
        self.k_tournament = k_tournament

        self.crossover_type = crossover_type
        self.crossover_prob = crossover_prob

        self.mutation_prob = mutation_prob

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
