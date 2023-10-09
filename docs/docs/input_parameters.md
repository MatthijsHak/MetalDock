# Input  

MetalDock is configured using a **'.ini'** file in which various parameters are set. This section provides an overview of the input parameters, categorized into the different headers of the **'.ini'** file. For examples, see **'input_examples'** of the GitHub repository.

## DEFAULT keywords 

---

### method
Describes the method which method to use in MetalDock.

**Default input:**  
method = dock  

**Valid values:**   
`dock` - Perform docking protocol for a single metal-organic compound.  
`mc`   - Perform Monte Carlo optimisation scheme   

---

### metal_symbol
The metal symbol of the metal atom in the metal-organic compound that you want to dock. The list of metals that er given as valid values are currently the only metals that MetalDock can dock.


**Default input:**  
metal_symbol = Ru

**Valid values:**  
`V`  - Vanadium is used as metal in MetalDock.  
`Cr` - Chromium is used as metal in MetalDock.   
`Co` - Cobalt is used as metal in MetalDock.   
`Ni` - Nickel is used as metal in MetalDock.    
`Cu` - Copper is used as metal in MetalDock.    
`Mo` - Molybdenum is used as metal in MetalDock.    
`Ru` - Ruthenium is used as metal in MetalDock.    
`Rh` - Rhodium is used as metal in MetalDock.    
`Pd` - Palladium is used as metal in MetalDock.    
`Re` - Rhenium is used as metal in MetalDock.   
`Os` - Osmium is used as metal in MetalDock.  
`Pt` - Platinum is used as metal in MetalDock.  

---

### parameter_file
The path to the parameter file used. If parameter file is specified as metal_dock.dat the internal parameters are used. If no parameter file is specified the internal parameters of MetalDock will be used. 


**Default input:**   
parameter_file = metal_dock.dat  

**Valid values:**    
`string`  - path to the parameter file 

---

### ncpu
The number of CPUs used in the quantum mechanical calculations.


**Default input:**  
ncpu = 1 

**Valid Values:**  
`integer` -  range[0,inf]  

---

## PROTEIN keywords 

---

### pdb_file
The path to the PDB file of the protein to which the metal-organic compound is docked. 

**Default input:**  
pdb_file = protein.pdb

**Valid values:**  
`string` - path to PDB file  

---

### pH
The pH of the system in which the docking experiment is performed.   

**Default input:**  
pH = 7.4

**Valid values:**  
`float` - range[0,14]

---

### clean_pdb
This keyword ensures that in the PDB file all atoms will be removed that are not from the protein. If a co-factor is present during docking set this keyword to False and remove all atoms yourself. 

**Default input:**  
clean_pdb = True

**Valid values:**  
`boolean` 

---

## METAL_COMPLEX Keywords 

---

### geom_opt 
This keyword activates a quantum mechanical geometry optimization of your metal-complex if set to True. 

**Default input:**   
geom_opt = True

**Valid values:**  
`boolean`

---

### xyz_file 
The path to the XYZ file of the metal-organic compound.

**Default input:**   
xyz_file = metal_complex.xyz

**Valid values:**  
`string` - path  

---

### charge
The total charge of the metal-organic compound. Spin is here defined as the number of spin-alpha electrons in excess of spin-beta electrons.

**Default input:**  
charge = 0 

**Valid values:**  
`integer`

---

### spin 
The total spin of the metal-organic compound. 

**Default input:**   
spin = 0

**Valid values:**  
`integer`  

---

### vacant_site
Specify whether there is a vancy in the first coordiantion sphere. For docking procedure where it is expected that the metal atom would coordinate directly to an atom of the residue, this keyword should be set to True. If there is no vanacy in the coordination sphere and the metal does not directly coordinated to the protein, this keyword should be set to False

**Default input:**  
vacant_site = True 

**Valid values:**  
`boolean`

---

## QM keywords 

---

### engine
The quantum mechanical engine that will be used for the single point calculations and geometry optimizations.

**Default input:**  
engine = ADF

**Valid values:**  
`ADF` -  the ADF engine will be used  
`Gaussian` - the Gaussian engine will be used   
`ORCA` - the ORCA engine will be used  

---

### basis_set (only ADF & Gaussian)
The basis set used for the DFT calculations.

**Default input:**  
basis_set = TZP

**Valid values:**  
`ADF` valid inputs can be found [here](https://www.scm.com/doc/ADF/Input/Basis_sets_and_atomic_fragments.html)   
`Gaussian` valid inputs can be found [here](https://gaussian.com/basissets/)  

---

### functional_type (only ADF & Gaussian)
The functional type used for the DFT calculations.

**Default input:**  
functional_type = GGA

**Valid values:**  
`ADF` valid inputs can be found [here](https://www.scm.com/doc/ADF/Input/Density_Functional.html)   
`Gaussian` valid inputs can be found [here](https://gaussian.com/dft/)  

---

### functional (only ADF & Gaussian)
The basis set used for the DFT calculations.

**Default input:**  
functional = PBE

**Valid values:**  
`ADF` valid inputs can be found [here](https://www.scm.com/doc/ADF/Input/Density_Functional.html)   
`Gaussian` valid inputs can be found [here](https://gaussian.com/dft/)  

---

### dispsersion (only ADF & Gaussian)
The basis set used for the DFT calculations.

**Default input:**  
dispsersion = Grimme3 BJDAMP

**Valid values:**  
`ADF` valid inputs can be found [here](https://www.scm.com/doc/ADF/Input/Density_Functional.html#dispersion-corrections)   
`Gaussian` valid inputs can be found [here](https://gaussian.com/dft/)  

---

### solvent (only ADF & Gaussian)
The basis set used for the DFT calculations.

**Default input:**  
solvent = water

**Valid values:**  
`ADF` valid inputs can be found [here](https://www.scm.com/doc/ADF/Input/COSMO.html)   
`Gaussian` valid inputs can be found [here](https://gaussian.com/basissets/)  

---

### orcasimpleinput (only ORCA)
This keyword accepts a string for the values where you in an ORCA run script would put a ! in front

**Default input:**  
orcasimpleinput = PBE def2-TZVP

**Valid values:**  
`ORCA` valid inputs can be found [here](https://sites.google.com/site/orcainputlibrary/home) 

---

### orcablocks (only ORCA)
This keywords accepts a string for the values where you in an ORCA run script would put a % in front. You should **NOT** specify here the number of CPUs, as that is already taken of with ncpu keyword.

**Default input:**  
orcablocks = 

**Valid values:**  
`ORCA` valid inputs can be found [here](https://sites.google.com/site/orcainputlibrary/home) 

---

## DOCKING keywords

### rmsd
Calculate the root mean squared deviation (RMSD) of the poses obtained from the docking procedure and the XYZ coordinates of the XYZ file used as input.

**Default input:**  
rmsd = False

**Valid values:**  
`boolean`

---

### dock_x
The x coordinate of the centre of the box used in the docking procedure.

**Default input:**  
dock_x = 0.0

**Valid values:**  
`float`

---

### dock_y
The y coordinate of the centre of the box used in the docking procedure.

**Default input:**  
dock_y = 0.0

**Valid values:**  
`float`

---

### dock_z
The x coordinate of the centre of the box used in the docking procedure.

**Default input:**  
dock_z = 0.0

**Valid values:**  
`float`

---

### box_size
This keyword specifies the box size side in Ångstrom of the box used in the docking procedure

**Default input:**  
box_size = 20

**Valid values:**  
`float`

---

### box_size
This keyword specifies the box size side in Ångstrom of the box used in the docking procedure

**Default input:**  
box_size = 20

**Valid values:**  
`float`

---

### scale_factor
This keyword scales the box to the volume of the metal-organic compound. 

**Default input:**  
box_size = 

**Valid values:**  
`float`

---

### random_pos
Randomize the initial positions of the atoms within the box before the docking procedure. 

**Default input:**  
box_size = True

**Valid values:**  
`boolean`

---

### ini_parameters
Keyword that needs to be activated if the well-depth parameters of the metal protein interaction are specified in the **'.ini'** file

**Default input:**  
ini_parameters = False

**Valid values:**  
`boolean`

---

### e_NA
The well-depth parameter of the Lennard-Jones interaction between the metal and the nitrogen hydrogen bond accepting atom type of the protein in kcal/mol. This function is only triggered when the standard keyword is set to False

**Default input:**  
box_size = 5.0

**Valid values:**  
`float`

---

### e_OA
The well-depth parameter of the Lennard-Jones interaction between the metal and the oxygen hydrogen bond accepting atom type of the protein in kcal/mol. This function is only triggered when the standard keyword is set to False

**Default input:**  
box_size = 5.0

**Valid values:**  
`float`

---

### e_SA
The well-depth parameter of the Lennard-Jones interaction between the metal and the sulphur hydrogen bond accepting atom type of the protein in kcal/mol. This function is only triggered when the standard keyword is set to False

**Default input:**  
box_size = 5.0

**Valid values:**  
`float`

---

### e_HD
The well-depth parameter of the Lennard-Jones interaction between the metal and the hydrogen bond donating atom type of the protein in kcal/mol. This function is only triggered when the standard keyword is set to False

**Default input:**  
box_size = 5.0

**Valid values:**  
`float`

---

### num_poses
Number of poses that are docked. 

**Default input:**  
num_poses = 10

**Valid values:**  
`int`

---

### ga_dock
If set to True the docking procedure will be performed with a genetic-algorithm.

**Default input:**  
ga_dock = True

**Valid values:**  
`boolean`

---

### ga_dock_pop_size
The number of individuals in the population of the genetic algorithm.

**Default input:**  
ga_dock_pop_size = 150

**Valid values:**  
`integer`

---

### ga_dock_num_evals
The maximum number of energy evluations performed during each genetic algorithm run. 

**Default input:**  
ga_dock_num_evals = 2500000

**Valid values:**  
`integer`

---

### ga_dock_num_generations
The maximum number of generations simulated during each genetic algorithm run.

**Default input:**  
ga_dock_num_generations = 27000

**Valid values:**  
`integer`

---

### ga_dock_elitism
The number of top individuals that are guaranteed to survive into the next generation.

**Default input:**  
ga_dock_elitism = 1

**Valid values:**  
`integer`

---

### ga_dock_mutation_rate
The probability that a particulare gene is mutated in the docking procedure of the genetic algorithm.

**Default input:**  
ga_dock_mutation_rate = 0.02

**Valid values:**  
`float` - [0,1]

---

### ga_dock_crossover_rate
Crossover rate is the expected number of pairs in the population that will exchange genetic material. Setting this value to 0 turns the GA into the evolutionary programming (EP) method, but EP would probably require a concomitant increase in the ga_mutation_rate in order to be effective.

**Default input:**  
ga_dock_crossover_rate = 0.80

**Valid values:**  
`float` - [0,1]

---

### ga_dock_window_size
The number of preceding generations to take into consideration when deciding the threshold for the worst individual in the current population.

**Default input:**  
ga_dock_window_size = 10

**Valid values:**  
`integer`

---

### sa_dock
If set to True the docking procedure will be performed with a simulated eannealing algorithm.

**Default input:**  
sa_dock = False

**Valid values:**  
`boolean`

---

### temp_reduction_factor
At the end of each cycle, the annealing temperature is multiplied by this factor, to give that of the next cycle. This must be positive but < 1 in order to cool the system. Gradual cooling is recommended, so as to avoid “simulated quenching”, which tends to trap systems into local minima.

**Default input:**  
temp_reduction_factor = 0.90

**Valid values:**  
`float`

---

### number_of_runs
Number of docking runs

**Default input:**  
number_of_runs = 50

**Valid values:**  
`integer`

---

### max_cycles
Number of temperature reduction cycles. 

**Default input:**  
max_cycles = 50

**Valid values:**  
`integer`

---

## MC keywords

### mc_steps
The number of steps taken in the Monte Carlo optimisation scheme

**Default input:**  
mc_steps = 250

**Valid values:**  
`integer`
