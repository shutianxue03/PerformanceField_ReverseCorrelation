#!/bin/bash
# ============================================================
# NOM pipeline launcher
# Author: Shutian Xue
# Last modified: 2026-01-08 (adapted for NYU Torch)
#
# Description:
#   Use flag_stage to choose which step to submit:
#     flag_stage=1 → run NOM_compIV (must be done BEFORE fitNOM)
#     flag_stage=2 → run NOM_fitNOM (requires outputs from compIV)
# ============================================================

set -euo pipefail

flag_stage=1    # 1 = run xx_compIV (MUST be done before running xx_fitNOM); 2 = run xx_fitNOM
nIter=100  # number of iterations
nJob=10 # number of jobs

for iJob in {1..nJob}; do # 
for isubj in {1..15}; do        # subject indices
  for iLocComb in {1..5}; do    # location combinations
    for iModelA in {1..2}; do      # see SX_RC1_setting (A1: core, A2: IO template)

      if [ "$flag_stage" -eq 1 ]; then
        echo "Submitting compIV: subj=${isubj}, loc=${iLocComb}, modelA=${iModelA}, nIter=${nIter}, iJob=${iJob}/${nJob}"

        sbatch \
          --job-name="compIV_S${isubj}_L${iLocComb}_A${iModelA}_Job${iJob}" \
          shell_NOM_compIV.sh "${isubj}" "${iLocComb}" "${iModelA}" "${nIter}" "${iJob}" "${nJob}" 

      elif [ "$flag_stage" -eq 2 ]; then
        for iModelB in {1..4}; do   # see SX_RC1_setting "namesModelB"
          echo "Submitting fitNOM: subj=${isubj}, loc=${iLocComb}, A=${iModelA}, B=${iModelB}, nIter=${nIter}, iJob=${iJob}/${nJob}"

          sbatch \
            --job-name="fitNOM_S${isubj}_L${iLocComb}_A${iModelA}_B${iModelB}_Job${iJob}" \
            shell_NOM_fitNOM.sh "${isubj}" "${iLocComb}" "${iModelA}" "${iModelB}" "${nIter}" "${iJob}" "${nJob}" 
        done
      else
        echo "ERROR: flag_stage must be 1 or 2 (current: ${flag_stage})" >&2
        exit 1
      fi

    done
  done
done
done
