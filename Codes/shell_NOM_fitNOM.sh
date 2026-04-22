#!/bin/bash
# ============================================================
# SLURM wrapper: OOD_NOM_Trialwise_fitNOM
# Author: Shutian Xue
# Last modified: 2026-01-08 (adapted for NYU Torch)
# 
# This script:
# 1. Runs a small performance evaluation script (optional monitoring)
# 2. Launches MATLAB and calls:
# loops over iModelB_fit=1:4, calling:
# OOD_NOM_Trialwise_fitNOM(nBasisORI, nBasisSF, isubj, iLocComb, flag_incluCrit, C_contribution, iModelA_fit, iModelB_fit, nIter, nJob, iJob)
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
echo "Subj=$1 Loc=$2 ModelA=$3 ModelB=1:4 nIter=$4 iJob=$6/$5"
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
iModelA_fit = $3;
nIter = $4;
nJob = $5;
iJob = $6;

nBasisORI = 6; % hard-coded, not used in this function
nBasisSF = 6; % hard-coded, not used in this function
flag_incluCrit = 1; % hard-coded; used; 1=include criterion as a free parameter
C_contribution = 0.5; % hard-coded; used; contribution of criterion to the model

for iModelB_fit = 1:4
    fprintf('\n--- ModelB=%d ---\n', iModelB_fit);
    OOD_NOM_Trialwise_fitNOM(nBasisORI, nBasisSF, isubj, iLocComb, flag_incluCrit, C_contribution, iModelA_fit, iModelB_fit, nIter, nJob, iJob)
end

disp('====================================');
disp(' OOD_NOM_Trialwise_fitNOM finished (all ModelB).');
disp('====================================');

exit
EOF

# Clean up temporary MATLAB preferences
rm -rf "$MATLAB_PREFDIR"

exit 0