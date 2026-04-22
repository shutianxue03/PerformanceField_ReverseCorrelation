#!/bin/bash
# ============================================================
# Bash/Slurm batch launcher for NOM pipeline
#
# Stage 1 (flag_stage=1): calls shell_NOM_compIV.sh
#   args: isubj iLocComb iModelA_fit nIter nJob iJob
#   hard-coded inside wrapper: nBasisORI=6, nBasisSF=6,
#     lambda_whiten=1, flag_regressType=2
#
# Stage 2 (flag_stage=2): calls shell_NOM_fitNOM.sh
#   args: isubj iLocComb iModelA_fit nIter nJob iJob
#   loops over iModelB_fit=1:4 internally
#   hard-coded inside wrapper: nBasisORI=6, nBasisSF=6,
#     flag_incluCrit=1, C_contribution=0.5
# ============================================================

set -euo pipefail

flag_stage=1
nIter=200
nJob=5

isubjList=($(seq 1 2))
iLocCombList=($(seq 1 2))
iModelA_fit_List=(1)
iModelB_fit_List=(1 3)  # looped inside shell_NOM_fitNOM.sh, not submitted separately

die(){ echo "ERROR: $*" >&2; exit 1; }

case "${flag_stage}" in
  1|2) ;;
  *) die "flag_stage must be 1 or 2 (current: ${flag_stage})" ;;
esac

for iJob in $(seq 1 "$nJob"); do
  for isubj in "${isubjList[@]}"; do
    for iLocComb in "${iLocCombList[@]}"; do
      for iModelA_fit in "${iModelA_fit_List[@]}"; do

        if [[ "${flag_stage}" -eq 1 ]]; then
          echo "compIV: job=${iJob}/${nJob}, subj=${isubj}, loc=${iLocComb}, A=${iModelA_fit}, nIter=${nIter}"

          sbatch \
            --job-name="compIV_S${isubj}_L${iLocComb}_A${iModelA_fit}_J${iJob}" \
            shell_NOM_compIV.sh "${isubj}" "${iLocComb}" "${iModelA_fit}" "${nIter}" "${nJob}" "${iJob}"

        else
          echo "fitNOM: job=${iJob}/${nJob}, subj=${isubj}, loc=${iLocComb}, A=${iModelA_fit}, B=1:4, nIter=${nIter}"

          sbatch \
            --job-name="fitNOM_S${isubj}_L${iLocComb}_A${iModelA_fit}_J${iJob}" \
            shell_NOM_fitNOM.sh "${isubj}" "${iLocComb}" "${iModelA_fit}" "${nIter}" "${nJob}" "${iJob}"
        fi

      done
    done
  done
done
