#!/bin/bash
#SBATCH --job-name=sbi_run
#SBATCH --time=08:00:00
#SBATCH --mem=32000
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=grossmannr@cbs.mpg.de
#SBATCH --output=/ptmp/grossmannr/somato_model/logs/sbi/%A/%x_%A.log
#SBATCH --error=/ptmp/grossmannr/somato_model/logs/sbi/%A/%x_%A.log

export PYTHONUNBUFFERED=TRUE

export SIMDIR=/ptmp/grossmannr/somato_model/results_grossmannr
export WDDIR=/u/grossmannr/SomatosensoryLaminarModel
# RESDIR must contain Figures/Main/eeg_results/.../roi_epochswise/group_roi_*_ses-elec_*.csv
# (the measured electrical-stimulation target data). Adjust if stored elsewhere.
export RESDIR=/ptmp/grossmannr/somato_model/results_grossmannr

srun /u/grossmannr/miniforge3/envs/pyrates__env/bin/python \
    /u/grossmannr/SomatosensoryLaminarModel/Simulations/run_sbi.py \
    --num-rounds 3 \
    --num-sims 3000 \
    --workers ${SLURM_CPUS_PER_TASK} \
    --seed 0
