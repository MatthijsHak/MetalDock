#!/bin/bash
#SBATCH -N 1
#SBATCH --tasks-per-node=32
#SBATCH -t 1:00:00
#SBATCH -p genoa 
#

module load 2023 
module load AMS/2023.101-intelmpi

/gpfs/home6/mhakkennes/.micromamba/bin/envs/atm8.1.2/bin/python -u /home/mhakkennes/MetalDock/metaldock -i input.ini -m dock
