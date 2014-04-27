#!/bin/bash

# SLURM submission script for python jobs

#SBATCH --nodes=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hgudjonson@uchicago.edu

# output will be nothing, just put into dummy file
#SBATCH --output=slurm.out

# default partition, 16 cores w/32 GB ram
#SBATCH --partition=sandyb

module load python

python $1
