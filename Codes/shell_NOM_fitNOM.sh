#!/bin/bash
# ============================================================
# SLURM wrapper: OOD_NOM_Trialwise_fitNOM
# Author: Shutian Xue
# Last modified: 2026-01-08 (adapted for NYU Torch)
# 
# This script:
# 1. Runs a small performance evaluation script (optional monitoring)
# 2. Launches MATLAB and calls:
# OOD_NOM_Trialwise_fitNOM(isubj, iLocComb, iModelA, iModelB, nIter)
# ============================================================

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:59:59
#SBATCH --output=zzz_NOM2_%j.out
#SBATCH --account=torch_pr_503_general

set -euo pipefail

#############################
# MATLAB setup
#############################
module purge
module load matlab/2025b

# Use a temporary MATLAB preference directory to avoid conflicts
export MATLAB_PREFDIR
export SLURM_ACCOUNT="torch_pr_503_general"
export SBATCH_ACCOUNT=${SLURM_ACCOUNT}
export SALLOC_ACCOUNT=${SLURM_ACCOUNT}
MATLAB_PREFDIR=$(mktemp -d -t matlab-XXXX)

echo
echo "SLURM job ID : $SLURM_JOB_ID"
echo "SLURM job name: $SLURM_JOB_NAME"
echo "Subj=$1 Loc=$2 ModelA=$3 ModelB=$4 nIter=$5 iJob=$7/$6"
echo

#############################
# Run MATLAB non-interactively
#############################
cat <<EOF | srun matlab -nodisplay -nosplash -nodesktop

disp('====================================');
disp(' Starting OOD_NOM_Trialwise_fitNOM...');
disp('====================================');

isubj = $1;
iLocComb = $2;
iModelA = $3;
iModelB=$4;
nIter = $5;
nJob = $6;
iJob = $7;

OOD_NOM_Trialwise_fitNOM(isubj, iLocComb, iModelA, iModelB, nIter, nJob, iJob);

disp('====================================');
disp(' OOD_NOM_Trialwise_fitNOM finished.');
disp('====================================');

exit
EOF

# Clean up temporary MATLAB preferences
rm -rf "$MATLAB_PREFDIR"

exit 0