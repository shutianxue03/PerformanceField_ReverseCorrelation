#!/bin/bash
# Created by Shutian Xue on 11/23/2025
# Last modified by Shutian Xue on 11/23/2025

#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=zzz_NOM_%j.out
#SBATCH --mail-user=vivanxuest@gmail.com
#SBATCH --mail-type=END

#########performance evaluation########
wk_dir=$(pwd)
bash $wk_dir/evaluation-performance/evaluation-performance.sh $wk_dir/evaluation-performance/
############end########################

module purge
module load matlab/2024b
 export MATLAB_PREFDIR=$(mktemp -d -t matlab-XXXX)

echo
echo "job name: $SLURM_JOB_NAME"

cat<<EOF | srun matlab -nodisplay

%===============
isubj=$1;
iLocComb=$2;
iModelA=$3;
nIterations=$4;

OOD_NOM_Trialwise(isubj, iLocComb, iModelA, nIterations)

%===============
exit
EOF
 
rm -rf $MATLAB_PREFDIR
exit
