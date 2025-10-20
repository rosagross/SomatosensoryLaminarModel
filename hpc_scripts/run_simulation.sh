#!/bin/bash
#SBATCH --job-name=prepro_step004
#SBATCH --time=00:30:00
#SBATCH --mem=8000
#SBATCH --cpus-per-task=5
#SBATCH --array=13-52
#SBATCH --mail-user=grossmannr@cbs.mpg.de
#SBATCH --output=/ptmp/grossmannr/nt_detect/logs/prepro/%A/%x_%A_%a.log
#SBATCH --error=/ptmp/grossmannr/nt_detect/logs/prepro/%A/%x_%A_%a.log

export PYTHONUNBUFFERED=TRUE

export RAWDIR=/ptmp/grossmannr/nt_detect/rawdir                    # Set environment variables
export DATADIR=/ptmp/grossmannr/nt_detect/datadir 
export RESDIR=/ptmp/grossmannr/nt_detect/results_grossmannr
export WDDIR=/u/grossmannr/NearThresholdDetection

srun /u/grossmannr/miniforge3/envs/pyrates_env/bin/python /u/grossmannr/SomatosensoryLaminarModel/simulation_main.py --coupling $SLURM_ARRAY_TASK_ID