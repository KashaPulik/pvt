#!/bin/bash

#SBATCH --nodes=1 --ntasks-per-node=8
#SBATCH --job-name=omptask
#SBATCH --partition=release

cd $SLURM_SUBMIT_DIR

export OMP_NUM_THREADS=8

./app