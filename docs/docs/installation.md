# Installation

To install MetalDock correctly, please follow the instructions below:

## Step 1: Create Conda Environment

First, make sure you have Miniconda or Anaconda installed. If not, you can find installation instructions at [Miniconda](https://docs.conda.io/en/latest/miniconda.html).

Once you have Conda installed, create a new Conda environment named "MetalDock" and use Python 3.8:

```bash
conda create -n MetalDock python=3.8
```

Activate the "MetalDock" environment:

```bash
conda activate MetalDock
```

## Step 2: Install Required Packages 

Next, install several required packages using Conda and Pip:

```bash
conda install -c conda-forge openbabel
```

```bash
pip install ase rdkit-pypi pdb2pqr networkx numpy plams ase pandas
```

## Step 3: Clone MetalDock repository 

Clone the MetalDock respository to your desired location.

``` bash
git clone https://github.com/MatthijsHak/MetalDock
```

## Step 4: Setup QM environment variables

MetalDock can run on three different quantum chemistry software packages:
### ADF 
To run ADF as QM engine the following environment variables need to be exported:
``` bash
export AMSHOME=/full/path/to/ams/ams2022
export AMSBIN=$AMSHOME/bin
export AMSRESOURCES=$AMSHOME/atomicdata
export SCMLICENSE=$AMSHOME/license.txt
```

Relativistic effects with ZORA. For more infromation see [SCM](https://www.scm.com/doc/ADF/index.html) site.

### Gaussian
To run Gaussian as QM engine the following environment variables need to be exported.
``` bash
export g16root=/full/path/to/gaussian/g16
```

Relativistic effects with DKH. For more infromation see [Gaussian](https://gaussian.com/) site

### ORCA (free)
To run ORCA as QM engine with correct parralelization the full path to the orca binary has to be exported with the following command:
``` bash
export ASE_ORCA_COMMAND='/full/path/to/orca/orca PREFIX.inp > PREFIX.out'
```

Relativistic effects can be selected with input. For more infromation see [ORCA](https://www.orcasoftware.de/tutorials_orca/index.html) site


## Step 5: Add MetalDock to PATH

If you want to use MetalDock conventiely from any directory, you can add its executable to your PATH. Replace /user/directory/MetalDock with the actual path to your MetalDock installation directory:

```bash
export PATH=$PATH:/user/defined/path/MetalDock
```

This completes the installation process for MetalDock. You can now run metaldock by:

```
metaldock -i input.ini 
```

