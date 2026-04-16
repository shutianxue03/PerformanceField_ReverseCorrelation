#!/usr/bin/env bash
# ============================================================
# Created by Shutian Xue on 07/16/2025
# Last modified on 04/09/2026
#
# Description:
# Submit simulation jobs across parameter combinations.
#
# OOD_sim signature:
#   OOD_sim(noiseCST, gaborCST, nTrials, ...
#           Nmul_true, Nadd_true, Nshared_true, ...
#           Cz_true, lambda_whiten, ...
#           flag_regressType, flag_incluCrit, C_contribution, ...
#           iModelB_sim, nIter, nBasisORI, nBasisSF)
# ============================================================

set -euo pipefail

for noiseCST in 0.2; do
  for gaborCST in $(seq 0.3 0.1 0.5); do
    for nTrials in 8000; do
      for Nmul_true in 0.3; do
        for Nadd_true in 0; do
          for Nshared_true in 0; do
            for Cz_true in -0.4; do
              for lambda_whiten in 1; do
                for flag_regressType in 2; do          # 1 = Univariate; 2 = Multivariate + smoothing
                  for flag_incluCrit in 1; do         # include criterion loss term or not
                    for C_contribution in 0; do       # weight of criterion loss
                      for iModelB_sim in 2; do        # 1 = Full; 2 = No Nshared; 3 = No Nmul; 4 = No Nadd
                        for nIter in 20; do
                          for nBasisORI in 6; do
                            for nBasisSF in 5; do

                                
                                job_name="gCST${gaborCST}_Nm${Nmul_true}_Na${Nadd_true}_Ns${Nshared_true}_Cz${Cz_true}_Bsim${iModelB_sim}"
    
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
                                  "${iModelB_sim}" \
                                  "${nIter}" \
                                  "${nBasisORI}" \
                                  "${nBasisSF}"
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