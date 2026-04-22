#!/bin/bash
# ============================================================
# SLURM wrapper: OOD_NOM_Trialwise_compIV
# Author: Shutian Xue
# Last modified: 2026-01-08 (adapted for NYU Torch)
#
# This script:
# 1. Runs a small performance evaluation script (optional monitoring)
# 2. Launches MATLAB and calls:
# OOD_NOM_Trialwise_compIV(nBasisORI, nBasisSF, isubj, iLocComb, lambda_whiten, flag_regressType, iModelA_fit, nIter, nJob, iJob)
# ============================================================

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=16G
#SBATCH --time=01:30:00
#SBATCH --output=zzz_NOM1_%j.out
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
echo "SLURM job ID : ${SLURM_JOB_ID:-LOCAL}"
echo "SLURM job name: ${SLURM_JOB_NAME:-shell_NOM_compIV}"
echo "Subj=$1 Loc=$2 ModelA=$3 nIter=$4 iJob=$6/$5"
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
iModelA_fit = $3;
nIter = $4;
nJob = $5;
iJob = $6;

nBasisORI = 6; # hard-coded, not used in this function
nBasisSF = 6; # hard-coded, not used in this function
lambda_whiten = 1; # hard-coded; used 
flag_regressType = 2; # hard-coded; used; 2=multivariate regression + smoothing

OOD_NOM_Trialwise_compIV(nBasisORI, nBasisSF, isubj, iLocComb, lambda_whiten, flag_regressType, iModelA_fit, nIter, nJob, iJob)

disp('====================================');
disp(' OOD_NOM_Trialwise_compIV finished.');
disp('====================================');

exit
EOF

# Clean up temporary MATLAB preferences
rm -rf "$MATLAB_PREFDIR"

exit 0