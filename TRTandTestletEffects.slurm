#!/bin/bash
#SBATCH -J TRT_mvm
#SBATCH -o jobname_%j.txt
#SBATCH -e jobname_%j.err
#SBATCH -o /N/project/Ar00157/Paper_3_MVM/Results/jobname_%j.txt
#SBATCH -e /N/project/Ar00157/Paper_3_MVM/Results/jobname_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --time=4-00:00:00
#SBATCH --mem=240G
#SBATCH --partition=general
#SBATCH -A general
#SBATCH --mail-user=mbvaldiv@iu.edu
#SBATCH --mail-type=BEGIN,FAIL,END

#Set up environment

module load r
export OMP_NUM_THREADS=128
export OPENBLAS_NUM_THREADS=128

#Change directories

cd /N/project/Ar00157/Paper_3_MVM
#Commands to run

srun --cpus-per-task 128 Rscript DataSimAnalysisPaper3.R

