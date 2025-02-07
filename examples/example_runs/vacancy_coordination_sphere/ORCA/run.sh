#!/bin/bash
#SBATCH -N 1
#SBATCH --tasks-per-node=32
#SBATCH -t 10:00:00
#SBATCH -p genoa 
#

module load 2024
module load ORCA/6.0.1-gompi-2024a-avx2

export ASE_ORCA_COMMAND='/sw/arch/RHEL9/EB_production/2024/software/ORCA/6.0.1-gompi-2024a-avx2/bin/orca PREFIX.inp > PREFIX.out'

/gpfs/home6/mhakkennes/.micromamba/bin/envs/atm8.1.2/bin/python -u /home/mhakkennes/MetalDock/metaldock -i input.ini -m dock
