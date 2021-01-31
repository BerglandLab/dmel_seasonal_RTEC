#!/usr/bin/env bash
#
#SBATCH -J runPredictionModelPerms # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0:20:00 ### 6 hours
#SBATCH --mem 10G
#SBATCH -o /scratch/aob2x/drosRTEC/slurmOutput/runPredictionModelPerms.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/drosRTEC/slurmOutput/runPredictionModelPerms.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

## run as:

module load intel/18.0 intelmpi/18.0 R/3.6.3

wd=/scratch/aob2x/drosRTEC/drosRTEC_revisions

#SLURM_ARRAY_TASK_ID=1

#Rscript ${wd}/predictability_model/perm.thermalLimitModel.univariate.R ${SLURM_ARRAY_TASK_ID}
#Rscript ${wd}/predictability_model/perm.thermalLimitModel.bivariate.R ${SLURM_ARRAY_TASK_ID}
Rscript ${wd}/predictability_model/perm.altModels.R ${SLURM_ARRAY_TASK_ID}
