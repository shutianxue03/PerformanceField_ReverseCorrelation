#!/bin/bash
# ============================================================
# SLURM wrapper: OOD_NOM_Trialwise_fitNOM
# Author: Shutian Xue
# Last modified: 2025-11-23
#
# Usage (from launcher script):
#   sbatch shell_NOM_fitNOM.sh isubj iLocComb iModelA iModelB nIterations
#
# This script:
#   1. Runs a small performance evaluation script (optional monitoring)
#   2. Launches MATLAB and calls:
#        OOD_NOM_Trialwise_fitNOM(isubj, iLocComb, iModelA, iModelB, nIterations)
# ============================================================

#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=zzz_NOM2_fitNOM_%j.out
#SBATCH --mail-user=vivanxuest@gmail.com
#SBATCH --mail-type=END

set -euo pipefail

#############################
# Performance evaluation hook
#############################
wk_dir=$(pwd)
bash "$wk_dir/evaluation-performance/evaluation-performance.sh" "$wk_dir/evaluation-performance/"

#############################
# MATLAB setup
#############################
module purge
module load matlab/2024b

# Use a temporary MATLAB preference directory to avoid conflicts
export MATLAB_PREFDIR
MATLAB_PREFDIR=$(mktemp -d -t matlab-XXXX)

echo
echo "SLURM job ID  : $SLURM_JOB_ID"
echo "SLURM job name: $SLURM_JOB_NAME"
echo "Subj=$1  Loc=$2  ModelA=$3  nIter=$5"
echo


#############################
# Run MATLAB non-interactively
#############################
cat <<EOF | srun matlab -nodisplay -nosplash -nodesktop

disp('====================================');
disp('  Starting OOD_NOM_Trialwise_fitNOM...');
disp('====================================');

isubj     = $1;
iLocComb  = $2;
iModelA   = $3;
iModelB=$4;
nIterations = $5;

OOD_NOM_Trialwise_fitNOM(isubj, iLocComb, iModelA, iModelB, nIterations);

disp('====================================');
disp('  OOD_NOM_Trialwise_fitNOM finished.');
disp('====================================');

exit
EOF

# Clean up temporary MATLAB preferences
rm -rf "$MATLAB_PREFDIR"

exit 0