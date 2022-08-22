class lig_par_dock:

    def __init__(self, parameter_file, box_size,  metal_symbol, metal_cap, r_OA, e_OA, r_SA, e_SA, r_HD, e_HD, r_NA, e_NA, r_N, e_N, r_M, e_M):

        self.parameter_file = parameter_file
        self.box_size = box_size

        self.metal_symbol = metal_symbol
        self.metal_cap = metal_cap

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
