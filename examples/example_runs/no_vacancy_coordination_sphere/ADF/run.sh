#!/bin/bash
#SBATCH -N 1
#SBATCH --tasks-per-node=32
#SBATCH -t 1:00:00
#SBATCH -p genoa 
#

module load 2023 
module load AMS/2023.101-intelmpi

~/.conda/envs/MetalDock/bin/python -u /home/MetalDock/metaldock -i input.ini
