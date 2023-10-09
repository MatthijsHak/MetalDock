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
conda install -c bioconda mgltools
```

```bash
pip install ase rdkit-pypi pdb2pqr networkx numpy plams ase
```

## Step 3: Clone MetalDock repository 

Clone the MetalDock respository to your desired location.

``` bash
git clone https://github.com/MatthijsHak/MetalDock
```

## Step 4: Add MetalDock to PATH

If you want to use MetalDock conventiely from any directory, you can add its executable to your PATH. Replace /user/directory/MetalDock with the actual path to your MetalDock installation directory:

```bash
export PATH=$PATH:/Users/matthijsh/Desktop/MetalDock
```

This completes the installation process for MetalDock. You can now run metaldock by:

```
metaldock -i input.ini 
```