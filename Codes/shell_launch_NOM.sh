#!/bin/bash
# ============================================================
# NOM pipeline launcher (NYU Torch / Slurm)
# Author: Shutian Xue
# Last modified: 2026-01-08
#
# Usage:
#   Set flag_stage to choose which step to submit:
#     flag_stage=1  -> submit NOM_compIV   (must run before fitNOM)
#     flag_stage=2  -> submit NOM_fitNOM   (requires compIV outputs)
# ============================================================

set -euo pipefail

flag_stage=1    # 1 = compIV; 2 = fitNOM
nIter=2         # iterations per job
nJob=2          # number of jobs

isubjList=($(seq 1 15)) 
iLocCombList=($(seq 1 5)) 
iModelAList=(1 2)
iModelBList=(1 2 3 4)


# ---- Helpers ----
die() { echo "ERROR: $*" >&2; exit 1; }

case "${flag_stage}" in
  1|2) ;;
  *) die "flag_stage must be 1 or 2 (current: ${flag_stage})" ;;
esac

# ---- Submit jobs ----
for ((iJob=1; iJob<=nJob; iJob++)); do
  for isubj in "${isubjList[@]}"; do
    for iLocComb in "${iLocCombList[@]}"; do
      for iModelA in "${iModelAList[@]}"; do

        if [[ "${flag_stage}" -eq 1 ]]; then
          echo "Submitting compIV: subj=${isubj}, loc=${iLocComb}, A=${iModelA}, nIter=${nIter}, job=${iJob}/${nJob}"

          sbatch \
            --job-name="compIV_S${isubj}_L${iLocComb}_A${iModelA}_J${iJob}" \
            shell_NOM_compIV.sh "${isubj}" "${iLocComb}" "${iModelA}" "${nIter}" "${nJob}" "${iJob}"

        elif [[ "${flag_stage}" -eq 2 ]]; then
          for iModelB in "${iModelBList[@]}"; do
            echo "Submitting fitNOM: subj=${isubj}, loc=${iLocComb}, A=${iModelA}, B=${iModelB}, nIter=${nIter}, job=${iJob}/${nJob}"

            sbatch \
              --job-name="fitNOM_S${isubj}_L${iLocComb}_A${iModelA}_B${iModelB}_J${iJob}" \
              shell_NOM_fitNOM.sh "${isubj}" "${iLocComb}" "${iModelA}" "${iModelB}" "${nIter}" "${nJob}" "${iJob}"
          done

        else
          echo "ERROR: flag_stage must be 1 or 2 (current: ${flag_stage})" >&2
          exit 1
        fi

      done
    done
  done
done
