#!/bin/bash

##################
# slurm settings #
##################

# where to put stdout / stderr
#SBATCH --output=/nfs/scratch01/abaud/htonnele/logs/%x/%A_%a.out
#SBATCH --error=/nfs/scratch01/abaud/htonnele/logs/%x/%A_%a.err

# time limit in minutes
#SBATCH --time=01:00:00

# queue # checkout with `sacctmgr show qos format=name%10,maxwall%14,maxtrespu%20,maxtres%20,maxjobspu%10`
#SBATCH --qos=vshort

# memory (MB)
#SBATCH --mem=10000

# job name
#SBATCH --job-name get_resGWAS

### job array directive
###SBATCH --array=10-43

# cpu slots
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14 # same as in R script


#################
# start message #
#################
echo [$(date +"%Y-%m-%d %H:%M:%S")] starting on $(hostname)


# make bash behave more robustly
set -e
set -u
set -o pipefail


###############################################
# slurm arrays are 0 based, sed line counting #
# starts from 1                               #
###############################################
#TASK_ID="$((SLURM_ARRAY_TASK_ID+1))"

###################
# set environment #
###################
# e.g.: module load Python/3.10.4-GCCcore-11.3.0

###############
# run command #
###############
cd /users/abaud/htonnele/git/lab/P50/16S/Figures
pwd 

echo "starting script in R"
Rscript porcupine_plot_dataPrep.R # your command
echo "done woth Rscript"

###############
# end message #
###############
echo [$(date +"%Y-%m-%d %H:%M:%S")] finished on $(hostname)

#sbatch run_porcupine_plot_dataPrep.sh