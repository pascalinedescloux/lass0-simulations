#!/bin/sh
#SBATCH --job-name=TV300
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --partition=shared
#SBATCH --time=06:00:00
	# time for R = 500:
	# smallGauss: < 1h
	# wideGauss: ~4:00
	# riboflavin: ~4:30
	# TV300: 
#SBATCH --mail-user=descloup
#SBATCH --mail-type=ALL
#SBATCH --array=0-10
	# smallGauss: 0-25?
	# wideGauss: 0-15
#SBATCH --export=type="TV300",ALL 

module load foss/2016b R/3.4.2
srun R CMD BATCH -${SLURM_ARRAY_TASK_ID} simuSparsityFixed.R simuRes.$type.sparsity.${SLURM_ARRAY_TASK_ID}.out