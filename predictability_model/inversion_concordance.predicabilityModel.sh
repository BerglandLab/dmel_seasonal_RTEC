#!/usr/bin/env bash
#
#SBATCH -J runPredictionModelPerms # A single job name for the array
#SBATCH --ntasks-per-node=5 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0:30:00 ### 6 hours
#SBATCH --mem 100G
#SBATCH -o /scratch/aob2x/drosRTEC/slurmOutput/runPredictionModelPerms.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/drosRTEC/slurmOutput/runPredictionModelPerms.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

## run as: sbatch --array=1-96 ${wd}/drosRTEC_revisions/predictability_model/inversion_concordance.predicabilityModel.sh
## sacct -j 19262543

## cat /scratch/aob2x/drosRTEC/slurmOutput/runPredictionModelPerms.19262442_3.err

module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3

wd=/scratch/aob2x/drosRTEC

#SLURM_ARRAY_TASK_ID=1
Rscript --vanilla ${wd}/drosRTEC_revisions/predictability_model/inversion_concordance.predicabilityModel.R ${SLURM_ARRAY_TASK_ID}
