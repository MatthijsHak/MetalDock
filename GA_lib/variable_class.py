class lig_par_dock:

    def __init__(self, parameter_file, box_size,  metal_symbol, metal_cap, num_generations, num_parents_mating, sol_per_pop, parent_selection_type, keep_parents, K_tournament, crossover_type, crossover_prob, mutation_type, mutation_prob, mutation_percent):

        self.parameter_file = parameter_file
        self.box_size = box_size

        self.metal_symbol = metal_symbol
        self.metal_cap = metal_cap

        self.sol_per_pop = sol_per_pop
        self.num_generations = num_generations
        self.num_parents_mating = num_parents_mating

        self.parent_selection_type = parent_selection_type
        self.keep_parents = keep_parents
        self.K_tournament = K_tournament

        self.crossover_type = crossover_type
        self.crossover_prob = crossover_prob

        self.mutation_type = mutation_type
        self.mutation_prob = mutation_prob
        self.mutation_percent = mutation_percent

