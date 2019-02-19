#!/bin/sh
#SBATCH --job-name=smallGauss
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --partition=shared
#SBATCH --time=02:00:00
	# time for R = 500:
	# smallGauss: < 1h
#SBATCH --mail-user=descloup
#SBATCH --mail-type=ALL
#SBATCH --array=0-20
#SBATCH --export=type="smallGauss",ALL 

module load foss/2016b R/3.4.2
srun R CMD BATCH -${SLURM_ARRAY_TASK_ID} simuSparsityFixed.R simuRes.$type.sparsity.${SLURM_ARRAY_TASK_ID}.out