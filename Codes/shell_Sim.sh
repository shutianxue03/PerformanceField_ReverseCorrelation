#!/usr/bin/env bash
# Created by Shutian Xue on 07/16/2025
# Last modified by Shutian Xue on 03/26/2026

#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --output=aaa_Sim_%j.out
#SBATCH --mail-user=vivanxuest@gmail.com
#SBATCH --mail-type=END

set -euo pipefail

module purge
module load matlab/2025b
export MATLAB_PREFDIR
MATLAB_PREFDIR="$(mktemp -d -t matlab-XXXX)"

wk_dir="$(pwd)"

if [ "$#" -ne 10 ]; then
    echo "Error: expected 10 input arguments, got $#"
    echo "Usage:"
    echo "  shell_Sim.sh noiseCST gaborCST nTrials Nmul_true Nadd_true Nshared_true Cz_true lambda_whiten iModelA_sim iModelB_sim"
    rm -rf "$MATLAB_PREFDIR"
    exit 1
fi

noiseCST="$1"
gaborCST="$2"
nTrials="$3"
Nmul_true="$4"
Nadd_true="$5"
Nshared_true="$6"
Cz_true="$7"
lambda_whiten="$8"
iModelA_sim="$9"
iModelB_sim="${10}"

cat <<EOF | srun matlab -nodisplay
disp('======================================');
disp('MATLAB job started');
fprintf('noiseCST      = %g\n', $noiseCST);
fprintf('gaborCST      = %g\n', $gaborCST);
fprintf('nTrials       = %d\n', $nTrials);
fprintf('Nmul_true     = %g\n', $Nmul_true);
fprintf('Nadd_true     = %g\n', $Nadd_true);
fprintf('Nshared_true  = %g\n', $Nshared_true);
fprintf('Cz_true       = %g\n', $Cz_true);
fprintf('lambda_whiten = %g\n', $lambda_whiten);
fprintf('iModelA_sim   = %d\n', $iModelA_sim);
fprintf('iModelB_sim   = %d\n', $iModelB_sim);
disp('======================================');

disp('Calling OOD_sim...');

noiseCST = $noiseCST;
gaborCST = $gaborCST;
nTrials = $nTrials;
Nmul_true = $Nmul_true;
Nadd_true = $Nadd_true;
Nshared_true = $Nshared_true;
Cz_true = $Cz_true;
lambda_whiten = $lambda_whiten;
iModelA_sim = $iModelA_sim;
iModelB_sim = $iModelB_sim;

OOD_sim(noiseCST, gaborCST, nTrials, Nmul_true, Nadd_true, Nshared_true, Cz_true, lambda_whiten, iModelA_sim, iModelB_sim);

disp('MATLAB job finished');
disp('======================================');
exit
EOF

echo "Cleaning up MATLAB_PREFDIR..."
rm -rf "$MATLAB_PREFDIR"

echo "SLURM job finished."
exit 0