#!/bin/bash
#SBATCH --job-name=bifurcation_sI026
#SBATCH --output=/data/p_02989/Modelling/grossmannr_wd/logs/%j_bifurcation.out
#SBATCH --error=/data/p_02989/Modelling/grossmannr_wd/logs/%j_bifurcation.err
#SBATCH --time=12:00:00           # Set a realistic wall time — job dies cleanly instead of hanging
#SBATCH --mem=16G                 # Adjust based on your model size
#SBATCH --cpus-per-task=4         # AUTO-07p is mostly serial, but PyRates compilation can use threads
#SBATCH --partition=standard # Replace with your cluster's partition name

# ── Environment ──────────────────────────────────────────────────────────────
export WDDIR="/data/p_02989/Modelling/mecozzi_wd/SomatosensoryLaminarModel"
export SIMDIR="/data/p_02989/Modelling/mecozzi_wd/Simulations"
export PYRATES_CACHE="1"

# Activate your conda environment
source /data/u_mecozzi_software/miniforge3/etc/profile.d/conda.sh
conda activate pyrates_project

# AUTO-07p needs this on PATH
export PATH="/data/u_mecozzi_software/miniforge3/envs/pyrates_project/auto-07p/bin:$PATH"

# ── Logging ──────────────────────────────────────────────────────────────────
echo "Job ID:     $SLURM_JOB_ID"
echo "Node:       $SLURM_NODELIST"
echo "Start time: $(date)"

# ── Run ──────────────────────────────────────────────────────────────────────
python main.py

echo "End time: $(date)"
