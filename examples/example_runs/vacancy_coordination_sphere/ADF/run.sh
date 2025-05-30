#!/bin/bash                                                                                                                                          
#SBATCH --nodes=1              # Number of nodes
#SBATCH --ntasks=16           # Total number of MPI processes (tasks)
#SBATCH --cpus-per-task=1      # Number of CPUs per ta
#SBATCH -t 20:00:00
#SBATCH -p genoa
#

module load 2023 
module load AMS/2023.101-intelmpi

/gpfs/home6/mhakkennes/.micromamba/bin/envs/atm8.1.2/bin/python -u /home/mhakkennes/MetalDock/metaldock -i input.ini -m dock
