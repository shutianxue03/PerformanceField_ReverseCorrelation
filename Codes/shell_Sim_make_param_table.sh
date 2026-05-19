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
printf "noiseCST\tgaborCST\tnTrials\tNmul_true\tNadd_true\tNshared_true\tcSDT_true\tlambda_whiten\tiModelB_sim\tnIter\tnBasisORI\tbasisWidthORI\tnBasisSF\tbasisWidthSF\n" > "${param_file}"

for gaborCST in 0.3; do
  for Nmul_true in 0.5; do
    for Nadd_true in 5; do
      for Nshared_true in 5; do
        for cSDT_true in 0; do
          for iModelB_sim in 6; do
            for nBasisORI in 4 5 6 7 8; do
              for basisWidthORI in .6 .8 1; do
                for nBasisSF in 4 5 6 7 8; do
                  for basisWidthSF in .4 .6 .8; do
                    for nIter in 100; do
                      for nTrials in 10000; do
                        for noiseCST in 0.2; do
                          for lambda_whiten in 0; do
                            printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
                              "${noiseCST}" \
                              "${gaborCST}" \
                              "${nTrials}" \
                              "${Nmul_true}" \
                              "${Nadd_true}" \
                              "${Nshared_true}" \
                              "${cSDT_true}" \
                              "${lambda_whiten}" \
                              "${iModelB_sim}" \
                              "${nIter}" \
                              "${nBasisORI}" \
                              "${basisWidthORI}" \
                              "${nBasisSF}" \
                              "${basisWidthSF}" >> "${param_file}"
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
