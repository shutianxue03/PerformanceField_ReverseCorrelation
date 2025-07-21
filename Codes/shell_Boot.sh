#!/bin/bash
# Created by Shutian Xue on 07/21/2025
# Last modified by Shutian Xue on 07/21/2025

#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=zzz_Boot_%j.out
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
# Parse input arguments from shell_all_boot.sh

OOD_boot_current($1, $2, $3, $4, $5, $6, $7, $8);

%===============
exit
EOF
 
rm -rf $MATLAB_PREFDIR
exit
