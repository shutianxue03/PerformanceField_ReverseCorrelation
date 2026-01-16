#!/bin/bash
# ============================================================
# SLURM wrapper: OOD_NOM_Trialwise_compIV
# Author: Shutian Xue
# Last modified: 2026-01-08 (adapted for NYU Torch)
#
# This script:
# 1. Runs a small performance evaluation script (optional monitoring)
# 2. Launches MATLAB and calls:
# OOD_NOM_Trialwise_compIV(isubj, iLocComb, iModelA, nIter)
# ============================================================

#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=zzz_NOM1_compIV_%j.out

set -euo pipefail

#############################
# MATLAB setup
#############################
module purge
module load matlab/2024b

# Use a temporary MATLAB preference directory to avoid conflicts
export MATLAB_PREFDIR
MATLAB_PREFDIR=$(mktemp -d -t matlab-XXXX)

echo
echo "SLURM job ID : $SLURM_JOB_ID"
echo "SLURM job name: $SLURM_JOB_NAME"
echo "Subj=$1 Loc=$2 ModelA=$3 nIter=$4 iJob=$5"
echo


#############################
# Run MATLAB non-interactively
#############################
cat <<EOF | srun matlab -nodisplay -nosplash -nodesktop

disp('====================================');
disp(' Starting OOD_NOM_Trialwise_compIV...');
disp('====================================');

isubj = $1;
iLocComb = $2;
iModelA = $3;
nIter = $4;
iJob = $5;
nJob = $6;

OOD_NOM_Trialwise_compIV(isubj, iLocComb, iModelA, nIter, iJob, nJob);

disp('====================================');
disp(' OOD_NOM_Trialwise_compIV finished.');
disp('====================================');

exit
EOF

# Clean up temporary MATLAB preferences
rm -rf "$MATLAB_PREFDIR"

exit 0