class lig_par_dock:

    def __init__(self, reference_docking, dock_x, dock_y, dock_z, rmsd, pdb_file_protein, name_protein, xyz_file_ligand, name_ligand, parameter_file, box_size,  metal_symbol, metal_cap, charge_ligand, spin_ligand,  basis_set, gga_functional, hybrid_functional, dispersion_correction, r_nitrogen_A, eps_nitrogen_A, r_nitrogen, eps_nitrogen):

        self.reference_docking = reference_docking

        self.dock_x = dock_x
        self.dock_y = dock_y
        self.dock_z = dock_z

        self.rmsd = rmsd

        self.pdb_file_protein = pdb_file_protein
        self.name_protein = name_protein

        self.xyz_file_ligand = xyz_file_ligand
        self.name_ligand = name_ligand

        self.parameter_file = parameter_file
        self.box_size = box_size

        self.metal_symbol = metal_symbol
        self.metal_cap = metal_cap
        self.charge_ligand = charge_ligand
        self.spin_ligand = spin_ligand

        self.basis_set = basis_set
        self.gga_functional = gga_functional
        self.hybrid_functional = hybrid_functional
        self.dispersion_correction = dispersion_correction

        self.r_nitrogen_A = r_nitrogen_A
        self.eps_nitrogen_A = eps_nitrogen_A

        self.r_nitrogen = r_nitrogen
        self.eps_nitrogen = eps_nitrogen

