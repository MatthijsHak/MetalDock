#!/bin/bash
#SBATCH -N 1
#SBATCH --tasks-per-node=32
#SBATCH -t 10:00:00
#SBATCH -p genoa 
#

module load 2022
module load ORCA/5.0.4-foss-2022a

export ASE_ORCA_COMMAND='/sw/arch/RHEL8/EB_production/2022/software/ORCA/5.0.4-foss-2022a/bin/orca PREFIX.inp > PREFIX.out'

~/.conda/envs/MetalDock/bin/python -u /home/MetalDock/metaldock -i input.ini
