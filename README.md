
# MetalDock 

MetalDock is an open software tool that can dock organometallic compounds to proteins with parameters that have been obtained with a genetic algorithm. The library consists of open docking software (AutoDock) and open quantum software (ORCA) and can be used by anyone. Three options are available:

1) Dock organometallic compounds
2) Train the genetic algorithm on a dataset
3) Test the parameters obtained from the genetic algorithm


## Installation





## Dock organometallic compounds
MetalDock has as main function to dock organometallic compounds easily and reproducible. The user only has to supply:
1) xyz file of the compound 
2) pdb file of the protein/DNA/biomolecule 
3) input file consisting of the various parameters (see example_files directory) 

The program will use QM software to optimize the geometry, if specified, and will extract CM5 charges via a single point calculation. The script can run on three different QM software packages:
# 1) ADF 
To run ADF as QM engine the following environment variables need to be exported:
``` bash
export AMSHOME=/full/path/to/ams/ams2022
export AMSBIN=$AMSHOME/bin
export AMSRESOURCES=$AMSHOME/atomicdata
export SCMLICENSE=$AMSHOME/license.txt
```

Relativistic effects with ZORA

# 2) Gaussian
To run Gaussian as QM engine the following environment variables need to be exported
``` bash
export g16root=/full/path/to/gaussian/g16
```

Relativistic effects with DKH

# 3) ORCA (free)
To run ORCA as QM engine with correct parralelization the full path to the orca binary has to be exported with the following command:
``` bash
export ASE_ORCA_COMMAND='/full/path/to/orca/orca PREFIX.inp > PREFIX.out'
```

Relativistic effects with ZORA

##
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


## Training procedure genetic algorithm 
The only requirements for an input are a xyz file of the compound, a pdb file of the protein and a input.ini file. 
```mermaid

graph TD;

%% Nodes

A([Start GA])

B([Initialize Parameter Sets])

C([Calculate Fitness Function])

D([Select Best Parameter Set])

E([Exchange of Parameters Between Best Sets])

F([Randomly Change Some Parameters])

G([New Parameter Sets Created])

H{{Exceeded Number of Generations? <br> OR <br> Fitness Function > Limit}}

I([Best Solution])

  

%% Links

A --> B

B --> C

C --> D

D --> E

E --> F

F --> G

G --> H

H -- NO --> C

H -- YES --> I

  

%% Class Definitions

classDef Class fill:#FFFFFF, stroke:#0045A9,stroke-width:2px, font-size:15px, font-family:palatino;

%% Assign Class Style to Nodes

class A,B,C,D,E,F,G,H,I Class;

  

%% Changing color of links [NOTE: Link arrows will remain black]

linkStyle default fill:none, font-family:palatino;
```
