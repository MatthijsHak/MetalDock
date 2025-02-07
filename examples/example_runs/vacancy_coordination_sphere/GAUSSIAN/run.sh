#!/bin/bash
#SBATCH -N 1
#SBATCH --tasks-per-node=32
#SBATCH -t 10:00:00
#SBATCH -p genoa 
#


module load 2024
module load Gaussian/g16.c02

/gpfs/home6/mhakkennes/.micromamba/bin/envs/atm8.1.2/bin/python -u /home/mhakkennes/MetalDock/metaldock -i input.ini -m dock
