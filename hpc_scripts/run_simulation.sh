#!/bin/bash
#SBATCH --job-name=simulation_run
#SBATCH --time=00:30:00
#SBATCH --mem=8000
#SBATCH --cpus-per-task=5
#SBATCH --array=10, 20, 30, 40, 50
#SBATCH --mail-user=grossmannr@cbs.mpg.de
#SBATCH --output=/ptmp/grossmannr/somato_model/logs/prepro/%A/%x_%A_%a.log
#SBATCH --error=/ptmp/grossmannr/somato_model/logs/prepro/%A/%x_%A_%a.log

export PYTHONUNBUFFERED=TRUE

export SIMDIR=/ptmp/grossmannr/somato_model/results_grossmannr
export WDDIR=/u/grossmannr/SomatosensoryLaminarModel

srun /u/grossmannr/miniforge3/envs/pyrates__env/bin/python /u/grossmannr/SomatosensoryLaminarModel/simulation_main.py --coupling $SLURM_ARRAY_TASK_ID