#!/usr/bin/env bash
# ============================================================
# Created by Shutian Xue on 07/16/2025
# Last modified on 04/21/2026
#
# Description:
# Submit simulation jobs across parameter combinations.
#
# OOD_sim signature:
#   OOD_sim(noiseCST, gaborCST, nTrials, ...
#           Nmul_true, Nadd_true, Nshared_true, ...
#           cSDT_true, lambda_whiten, ...
#           flag_regressType, flag_incluCrit, C_contribution, ...
#           iModelB_sim, nIter, nBasisORI, nBasisSF)
# ============================================================

set -euo pipefail

for gaborCST in 0.5; do
  for Nmul_true in 0.9; do
    for Nadd_true in 9; do
      for Nshared_true in 9; do
        for cSDT_true in 0; do
          for iModelB_sim in 1 2 3 4 5 6 7; do  # 1 = Full; 2 = No Nmul; 3 = No Nadd; 4 = No Nshared; 5 = Nmul-only; 6 = Nadd-only; 7 = Nshared-only
            for nIter in 20; do
              for nTrials in 10000; do
                for noiseCST in 0.2; do
                  for lambda_whiten in 0; do
                    for flag_regressType in 2; do  # 1 = Univariate; 2 = Multivariate + smoothing
                      for flag_incluCrit in 1; do  # include criterion loss term or not
                        for C_contribution in 0; do  # weight of criterion loss
                          for nBasisORI in 6; do
                            for nBasisSF in 5; do

                              job_name="gCST${gaborCST}_Nm${Nmul_true}_Na${Nadd_true}_Ns${Nshared_true}_Cz${cSDT_true}_Bsim${iModelB_sim}"

                              sbatch \
                                --job-name="${job_name}" \
                                ./shell_Sim.sh \
                                "${noiseCST}" \
                                "${gaborCST}" \
                                "${nTrials}" \
                                "${Nmul_true}" \
                                "${Nadd_true}" \
                                "${Nshared_true}" \
                                "${cSDT_true}" \
                                "${lambda_whiten}" \
                                "${flag_regressType}" \
                                "${flag_incluCrit}" \
                                "${C_contribution}" \
                                "${iModelB_sim}" \
                                "${nIter}" \
                                "${nBasisORI}" \
                                "${nBasisSF}"

                              # Print job name for tracking
                              echo "${job_name}"

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
done

echo "======================================"
echo "All jobs submitted."
echo "======================================"
