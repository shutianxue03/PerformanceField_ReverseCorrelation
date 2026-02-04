#!/bin/bash
# ============================================================
# MATLAB loop -> Bash/Slurm launcher
# Mirrors:
#   flag_stage, nIter, nJob
#   isubjList=1:15, iLocCombList=1:5, iModelAList=[1,2], iModelBList=[1..4]
#   for iJob=1:nJob ...
# ============================================================

set -euo pipefail

flag_stage=xx
nIter=200
nJob=5

isubjList=($(seq 1 15))
iLocCombList=($(seq 1 5))
iModelAList=(1 2)
iModelBList=(1 2 3 4)

die(){ echo "ERROR: $*" >&2; exit 1; }

case "${flag_stage}" in
  1|2) ;;
  *) die "flag_stage must be 1 or 2 (current: ${flag_stage})" ;;
esac

for iJob in $(seq 1 "$nJob"); do
  for isubj in "${isubjList[@]}"; do
    for iLocComb in "${iLocCombList[@]}"; do
      for iModelA in "${iModelAList[@]}"; do

        if [[ "${flag_stage}" -eq 1 ]]; then
          echo "compIV: job=${iJob}/${nJob}, subj=${isubj}, loc=${iLocComb}, A=${iModelA}, nIter=${nIter}"

          sbatch \
            --job-name="compIV_S${isubj}_L${iLocComb}_A${iModelA}_J${iJob}" \
            shell_NOM_compIV.sh "${isubj}" "${iLocComb}" "${iModelA}" "${nIter}" "${nJob}" "${iJob}"

        else
          for iModelB in "${iModelBList[@]}"; do
            echo "fitNOM: job=${iJob}/${nJob}, subj=${isubj}, loc=${iLocComb}, A=${iModelA}, B=${iModelB}, nIter=${nIter}"

            sbatch \
              --job-name="fitNOM_S${isubj}_L${iLocComb}_A${iModelA}_B${iModelB}_J${iJob}" \
              shell_NOM_fitNOM.sh "${isubj}" "${iLocComb}" "${iModelA}" "${iModelB}" "${nIter}" "${nJob}" "${iJob}"
          done
        fi

      done
    done
  done
done
