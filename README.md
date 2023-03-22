
# MetalDock 

MetalDock is an open software docking tool that can dock several organometallic compounds in an easy and reproducible manner to biomolecules, proteins or DNA.

## Citation

Please cite the following paper if you decide to use MetalDock in any project:


## Installation

To install MetalDock, follow these steps:

1. Download the latest release from the MetalDock GitHub repository.
2. Extract the files to a folder on your computer.
3. Install the required dependencies by running 'pip install .' from the command line
4. MetalDock is now ready to use and can be run by typing:
```python
python main.py -i input.ini 
```

## Usage
To dock organometallic compounds, the user only has to supply:
1) xyz file of the compound. 
2) pdb file of the protein/DNA/biomolecule.
3) input file consisting of the various parameters (see input_examples directory).

MetalDock can run on three different quantum chemistry software packages:
### 1) ADF 
To run ADF as QM engine the following environment variables need to be exported:
``` bash
export AMSHOME=/full/path/to/ams/ams2022
export AMSBIN=$AMSHOME/bin
export AMSRESOURCES=$AMSHOME/atomicdata
export SCMLICENSE=$AMSHOME/license.txt
```

Relativistic effects with ZORA.

### 2) Gaussian
To run Gaussian as QM engine the following environment variables need to be exported.
``` bash
export g16root=/full/path/to/gaussian/g16
```

Relativistic effects with DKH.

### 3) ORCA (free)
To run ORCA as QM engine with correct parralelization the full path to the orca binary has to be exported with the following command:
``` bash
export ASE_ORCA_COMMAND='/full/path/to/orca/orca PREFIX.inp > PREFIX.out'
```

Relativistic effects can be selected with input.

## Docking workflow
The program will protonate the protein/DNA/biomolecule at the specified pH. It will then generate two .pdbqt files that will be used in the AutoDock4.2 scheme. The parameter file consisting of the added metal parameters will be automatically generated, and our own derived parameters will be used. If desired, the standard parameters can be overwritten.

A workflow for the docking procedure is schematically given below.
```mermaid

graph TD;

%% Nodes

A([Organometallic Compound])

B([Geometry Optimization?])

C([Geometry Optimized?])

D([Single Point])

E([CM5 Charge])

G([Create .pdbqt])
  
H([Insert Charge in .pdbqt])
  
A1([Protein])

B1([Protonate at Correct pH])

C1([Clean Protein?])

D1([Delete Water and HETATM])

E1([Create .pdbqt])

A2([User Input Docking])

B2([Create .dat w/ Metal Parameters])

C2([Create .gpf & .dpf Files])

A3([Start Docking])

B3([Calculate RMSD?])

C3([Use Input as Reference])

D3([Write Docked Poses])

%% Links

A --> B

B --  Yes --> C

B --  No --> D

C -- Yes --> D

D --> E


E --> G

G --> H

H --> A3

A1 --> B1

B1 --> C1
  
C1 -- Yes --> D1

D1 --> E1

C1 -- No --> E1

E1 --> A3

A2 --> B2

B2 --> C2

C2 --> A3

A3 --> B3

B3 -- Yes --> C3

C3 --> D3

B3 -- No --> D3

%% Class Definitions

classDef Class fill:#FFFFFF, stroke:#0045A9,stroke-width:2px, font-size:15px, font-family:palatino;

%% Assign Class Style to Nodes

class A,B,C,D,E,F,G,H,I Class;

class A1,B1,C1,D1,E1 Class;

class A2,B2,C2 Class;

class A3,B3,C3,D3 Class;

%% Changing color of links [NOTE: Link arrows will remain black]

linkStyle default fill:none, font-family:palatino;
```


## Monte Carlo optimisation scheme 
