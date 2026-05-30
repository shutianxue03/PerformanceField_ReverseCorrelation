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
printf "noiseCST\tgaborCST\tnTrials\tNmul_true\tNadd_true\tNshared_true\tcSDT_true\tiModelB_sim\tnIter\n" > "${param_file}"

for gaborCST in 0.3; do
  for Nmul_true in $(seq 0.2 0.2 1.2); do
    for Nadd_true in $(seq 0.5 0.5 3); do
      for Nshared_true in $(seq 0.5 0.5 3); do
        for cSDT_true in -0.2 0 0.2; do
          for iModelB_sim in $(seq 1 1 7); do
            for nIter in 50; do
              for nTrials in 4000; do
                for noiseCST in 0.2; do
                  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
                    "${noiseCST}" \
                    "${gaborCST}" \
                    "${nTrials}" \
                    "${Nmul_true}" \
                    "${Nadd_true}" \
                    "${Nshared_true}" \
                    "${cSDT_true}" \
                    "${iModelB_sim}" \
                    "${nIter}" >> "${param_file}"
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
