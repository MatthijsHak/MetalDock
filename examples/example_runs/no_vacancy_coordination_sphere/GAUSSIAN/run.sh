#!/bin/bash
#SBATCH -N 1
#SBATCH --tasks-per-node=32
#SBATCH -t 10:00:00
#SBATCH -p genoa 
#

module load 2022 
module load Gaussian/g16.c02

~/.conda/envs/MetalDock/bin/python -u /home/MetalDock/metaldock -i input.ini
