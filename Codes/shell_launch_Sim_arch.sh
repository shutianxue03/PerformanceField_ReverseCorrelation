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

# Resolve paths robustly under Slurm.
# Note: when this launcher itself is submitted via sbatch, BASH_SOURCE may
# point to a temporary copy under /opt/slurm/... . Prefer SLURM_SUBMIT_DIR.
if [[ -n "${SLURM_SUBMIT_DIR:-}" && -f "${SLURM_SUBMIT_DIR}/shell_Sim.sh" ]]; then
  script_dir="${SLURM_SUBMIT_DIR}"
  project_root="$(cd "${script_dir}/.." && pwd)"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" && -f "${SLURM_SUBMIT_DIR}/Codes/shell_Sim.sh" ]]; then
  script_dir="${SLURM_SUBMIT_DIR}/Codes"
  project_root="${SLURM_SUBMIT_DIR}"
else
  script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  project_root="$(cd "${script_dir}/.." && pwd)"
fi

if [[ ! -f "${script_dir}/shell_Sim.sh" ]]; then
  echo "Error: cannot find ${script_dir}/shell_Sim.sh" >&2
  echo "Hint: submit from PF_RC root (sbatch Codes/shell_launch_Sim.sh) or from Codes (sbatch shell_launch_Sim.sh)." >&2
  exit 1
fi

for gaborCST in 0.3 0.4 0.5; do
  for Nmul_true in 0.2 0.4 0.6 0.8; do
    for Nadd_true in 2 4 6 8; do
      for Nshared_true in 2 4 6 8; do
        for cSDT_true in -0.4 -0.2 0 0.2; do
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
                                --chdir="${project_root}" \
                                "${script_dir}/shell_Sim.sh" \
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
