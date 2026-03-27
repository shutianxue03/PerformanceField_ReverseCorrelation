#!/usr/bin/env bash
# ============================================================
# Created by Shutian Xue on 07/16/2025
# Last modified on 03/26/2026
# Description: Submit simulation jobs for different models and parameters
# ============================================================

set -euo pipefail

for noiseCST in 0.2; do
  for gaborCST in 0.2; do
    for nTrials in 8000; do
      for Nmul_true in 0; do
        for Nadd_true in 0; do
          for Nshared_true in 0; do
            for Cz_true in 0; do
              for lambda_whiten in 0.5; do
                for iModelA_sim in 1; do
                  for iModelB_sim in 1; do

                    job_name="Sim_A${iModelA_sim}_B${iModelB_sim}_Nmul${Nmul_true}_Nadd${Nadd_true}_Nshared${Nshared_true}_Cz${Cz_true}_L${lambda_whiten}"
                    args="${noiseCST} ${gaborCST} ${nTrials} ${Nmul_true} ${Nadd_true} ${Nshared_true} ${Cz_true} ${lambda_whiten} ${iModelA_sim} ${iModelB_sim}"

                    echo "  job_name      = ${job_name}"
                    echo "  args          = ${args}"

                    sbatch \
                      --job-name="${job_name}" \
                      shell_Sim.sh \
                      "${noiseCST}" \
                      "${gaborCST}" \
                      "${nTrials}" \
                      "${Nmul_true}" \
                      "${Nadd_true}" \
                      "${Nshared_true}" \
                      "${Cz_true}" \
                      "${lambda_whiten}" \
                      "${iModelA_sim}" \
                      "${iModelB_sim}"

                    echo "Submitted ${job_name}"
                  done
                done
              done
            done
          done
        done
      done
    done
  done
done

echo "======================================"
echo "All jobs submitted."
echo "======================================"