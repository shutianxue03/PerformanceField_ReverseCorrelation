#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Create parameter table for OOD_NOM_Human array jobs.
# One line = one Human pipeline condition.
# ============================================================

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "${script_dir}/.." && pwd)"

param_dir="${project_root}/Data_job_params"
mkdir -p "${param_dir}"

param_file="${param_dir}/OOD_Human_params.tsv"
subject_ids=({1..15})
location_ids=(1 2 3 4 5)
norm_dv_flags=(1)
n_iters_all=(200)
n_jobs_all=(5)

# Header
printf "isubj\tiLocComb\tflag_normDV\tnIter\tnJob\tiJob\n" > "${param_file}"

for isubj in "${subject_ids[@]}"; do
  for iLocComb in "${location_ids[@]}"; do
    for flag_normDV in "${norm_dv_flags[@]}"; do
      for nIter in "${n_iters_all[@]}"; do
        for nJob in "${n_jobs_all[@]}"; do
          for ((iJob = 1; iJob <= nJob; iJob++)); do
            printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
              "${isubj}" \
              "${iLocComb}" \
              "${flag_normDV}" \
              "${nIter}" \
              "${nJob}" \
              "${iJob}" >> "${param_file}"
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
