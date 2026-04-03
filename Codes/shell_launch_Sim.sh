#!/usr/bin/env bash
# ============================================================
# Created by Shutian Xue on 07/16/2025
# Last modified on 03/30/2026
# Description: Submit simulation jobs for different models and parameters
# Updated to match OOD_sim signature:
# OOD_sim(noiseCST, gaborCST, nTrials, Nmul_true, Nadd_true, Nshared_true,
#         Cz_true, lambda_whiten, flag_regressType, flag_incluCrit,
#         C_contribution, iModelA_sim, iModelB_sim, nIter)
# ============================================================

set -euo pipefail

for noiseCST in 0.2; do
  for gaborCST in 0.2; do
    for nTrials in 8000; do
      for Nmul_true in 0; do
        for Nadd_true in 0; do
          for Nshared_true in 0; do
            for Cz_true in $(seq -0.5 0.1 0.2); do
              for lambda_whiten in $(seq 0 0.2 1); do
                for flag_regressType in 2; do            # 1=Univariate; 2=Multi+smoothing
                  for flag_incluCrit in 1; do          # include criterion loss term or not
                    for C_contribution in 0; do        # weight of criterion loss
                      for iModelA_sim in 1; do
                        for iModelB_sim in 1; do
                          for nIter in 20; do

                            job_name="Sim_A${iModelA_sim}_B${iModelB_sim}_Rg${flag_regressType}_Crit${flag_incluCrit}_Cw${C_contribution}_Nm${Nmul_true}_Na${Nadd_true}_Ns${Nshared_true}_Cz${Cz_true}_L${lambda_whiten}"

                            sbatch \
                              --job-name="${job_name}" \
                              ./shell_Sim.sh \
                              "${noiseCST}" \
                              "${gaborCST}" \
                              "${nTrials}" \
                              "${Nmul_true}" \
                              "${Nadd_true}" \
                              "${Nshared_true}" \
                              "${Cz_true}" \
                              "${lambda_whiten}" \
                              "${flag_regressType}" \
                              "${flag_incluCrit}" \
                              "${C_contribution}" \
                              "${iModelA_sim}" \
                              "${iModelB_sim}" \
                              "${nIter}"

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
      done
    done
  done
done

echo "======================================"
echo "All jobs submitted."
echo "======================================"
