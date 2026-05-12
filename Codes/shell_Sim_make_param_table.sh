#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Create parameter table for OOD_sim array jobs.
# One line = one simulation condition.
# ============================================================

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "${script_dir}/.." && pwd)"

param_dir="${project_root}/Data_job_params"
mkdir -p "${param_dir}"

param_file="${param_dir}/OOD_sim_params.tsv"

# Header
printf "noiseCST\tgaborCST\tnTrials\tNmul_true\tNadd_true\tNshared_true\tcSDT_true\tlambda_whiten\tflag_regressType\tflag_incluCrit\tC_contribution\tiModelB_sim\tnIter\tnBasisORI\tnBasisSF\n" > "${param_file}"

for gaborCST in 0.3; do
  for Nmul_true in 0.2 0.4 0.6 0.8; do
    for Nadd_true in 2 4 6 8; do
      for Nshared_true in 2 4 6 8; do
        for cSDT_true in -0.2 0 0.2; do
          for iModelB_sim in 1 2 3 4 5 6 7; do
            for nIter in 20; do
              for nTrials in 10000; do
                for noiseCST in 0.2; do
                  for lambda_whiten in 0; do
                    for flag_regressType in 2; do
                      for flag_incluCrit in 1; do
                        for C_contribution in 0; do
                          for nBasisORI in 6; do
                            for nBasisSF in 5; do
                              printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
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
                                "${nBasisSF}" >> "${param_file}"
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

n_jobs=$(( $(wc -l < "${param_file}") - 1 ))

echo "Wrote parameter table:"
echo "${param_file}"
echo "Number of jobs: ${n_jobs}"
echo
echo "Submit with:"
echo "sbatch --array=1-${n_jobs}%100 Codes/shell_array_Sim.sh"