#!/usr/bin/env bash
# ============================================================
# Created by Shutian Xue on 07/16/2025
# Last modified by Shutian Xue on 04/04/2026
# Description: Run one simulation job in Slurm
# Updated to match OOD_sim signature.
# ============================================================

#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=aaa_Sim_%j.out
#SBATCH --mail-user=vivanxuest@gmail.com
#SBATCH --mail-type=END

set -euo pipefail

module purge
module load matlab/2025b

export MATLAB_PREFDIR
MATLAB_PREFDIR="$(mktemp -d -t matlab-XXXX)"
trap 'rm -rf "$MATLAB_PREFDIR"' EXIT

wk_dir="$(pwd)"

noiseCST="$1"
gaborCST="$2"
nTrials="$3"
Nmul_true="$4"
Nadd_true="$5"
Nshared_true="$6"
Cz_true="$7"
lambda_whiten="$8"
flag_regressType="$9"
flag_incluCrit="${10}"
C_contribution="${11}"
iModelB_sim="${12}"
nIter="${13}"
nBasisORI="${14}"
nBasisSF="${15}"

cat <<EOF2 | srun matlab -nodisplay
try
    disp('======================================');
    disp('MATLAB job started');

    wk_dir = '$wk_dir';
    addpath(genpath(wk_dir));

    noiseCST = $noiseCST;
    gaborCST = $gaborCST;
    nTrials = $nTrials;
    Nmul_true = $Nmul_true;
    Nadd_true = $Nadd_true;
    Nshared_true = $Nshared_true;
    Cz_true = $Cz_true;
    lambda_whiten = $lambda_whiten;
    flag_regressType = $flag_regressType;
    flag_incluCrit = $flag_incluCrit;
    C_contribution = $C_contribution;
    iModelB_sim = $iModelB_sim;
    nIter = $nIter;
    nBasisORI = $nBasisORI;
    nBasisSF = $nBasisSF;

    fprintf('noiseCST         = %g\\n', noiseCST);
    fprintf('gaborCST         = %g\\n', gaborCST);
    fprintf('nTrials          = %g\\n', nTrials);
    fprintf('Nmul_true        = %g\\n', Nmul_true);
    fprintf('Nadd_true        = %g\\n', Nadd_true);
    fprintf('Nshared_true     = %g\\n', Nshared_true);
    fprintf('Cz_true          = %g\\n', Cz_true);
    fprintf('lambda_whiten    = %g\\n', lambda_whiten);
    fprintf('flag_regressType = %g\\n', flag_regressType);
    fprintf('flag_incluCrit   = %g\\n', flag_incluCrit);
    fprintf('C_contribution   = %g\\n', C_contribution);
    fprintf('iModelB_sim      = %g\\n', iModelB_sim);
    fprintf('nIter            = %g\\n', nIter);
    fprintf('nBasisORI        = %g\n', nBasisORI);
    fprintf('nBasisSF         = %g\n', nBasisSF);

    OOD_sim( ...
        noiseCST, ...
        gaborCST, ...
        nTrials, ...
        Nmul_true, ...
        Nadd_true, ...
        Nshared_true, ...
        Cz_true, ...
        lambda_whiten, ...
        flag_regressType, ...
        flag_incluCrit, ...
        C_contribution, ...
        iModelB_sim, ...
        nIter, ...
        nBasisORI, ...
        nBasisSF);

    disp('MATLAB job finished');
    disp('======================================');
catch ME
    disp(getReport(ME, 'extended'));
    exit(1);
end
exit
EOF2

echo "SLURM job finished."
exit 0
