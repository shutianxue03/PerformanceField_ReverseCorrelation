#!/bin/bash
# ============================================================
# NOM pipeline launcher
# Author: Shutian Xue
# Last modified: 2025-11-23
#
# Description:
#   Use flag_stage to choose which step to submit:
#     flag_stage=1 → run NOM_compIV (must be done BEFORE fitNOM)
#     flag_stage=2 → run NOM_fitNOM (requires outputs from compIV)
# ============================================================

set -euo pipefail

flag_stage=1    # 1 = run xx_compIV (MUST be done before running xx_fitNOM); 2 = run xx_fitNOM
nBoot=100  # number of bootstraps

for isubj in {1..15}; do        # subject indices
  for iLocComb in {1..8}; do    # location combinations
    for iModelA in {1..2}; do      # see SX_RC1_setting (A1: core, A2: IO template)

      if [ "$flag_stage" -eq 1 ]; then
        echo "Submitting compIV: subj=${isubj}, loc=${iLocComb}, modelA=${iModelA}, nIter=${nBoot}"

        sbatch \
          --job-name="compIV_S${isubj}_L${iLocComb}_A${iModelA}" \
          shell_NOM_compIV.sh "${isubj}" "${iLocComb}" "${iModelA}" "${nBoot}"

      elif [ "$flag_stage" -eq 2 ]; then
        for iModelB in {1..7}; do   # see SX_RC1_setting (B1–B7)
          echo "Submitting fitNOM: subj=${isubj}, loc=${iLocComb}, A=${iModelA}, B=${iModelB}, nIter=${nBoot}"

          sbatch \
            --job-name="fitNOM_S${isubj}_L${iLocComb}_A${iModelA}_B${iModelB}" \
            shell_NOM_fitNOM.sh "${isubj}" "${iLocComb}" "${iModelA}" "${iModelB}" "${nBoot}"
        done
      else
        echo "ERROR: flag_stage must be 1 or 2 (current: ${flag_stage})" >&2
        exit 1
      fi

    done
  done
done
