#!/bin/sh
#SBATCH --job-name=qutMC
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --partition=shared
#SBATCH --time=00:10:00
		# time for MCrep = 1000:
		# smallGauss: <10 minutes
		# wideGauss: <15 minutes
#SBATCH --mail-user=descloup
#SBATCH --mail-type=ALL
#SBATCH --array=1-10
#SBATCH --export=type="SmallGauss",ALL 


module load foss/2016b R/3.4.2
srun R CMD BATCH -${SLURM_ARRAY_TASK_ID} do.qut.MC.array.R do.qut.MC.$type.array${SLURM_ARRAY_TASK_ID}.out