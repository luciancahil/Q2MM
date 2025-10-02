#!/bin/sh

#SBATCH --time=10:00:00
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=6
#SBATCH --mem=128G
#SBATCH --job-name=Forces
#SBATCH --account=st-singha53-1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=royhe@student.ubc.ca
#SBATCH --output=forces.txt
#SBATCH --error=forces_error.txt


bash run_bo.sh open_mm_forces
