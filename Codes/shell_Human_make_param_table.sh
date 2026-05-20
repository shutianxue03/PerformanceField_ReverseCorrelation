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

# Header
printf "isubj\tiLocComb\tlambda_whiten\tnIter\n" > "${param_file}"

for isubj in $(seq 1 2 15); do
  for iLocComb in 1 2 3 4 5 6 7 8; do
    for lambda_whiten in 0; do
      for nIter in 1000; do
        printf "%s\t%s\t%s\t%s\n" \
          "${isubj}" \
          "${iLocComb}" \
          "${lambda_whiten}" \
          "${nIter}" >> "${param_file}"
      done
    done
  done
done

n_jobs=$(( $(wc -l < "${param_file}") - 1 ))

echo "Wrote parameter table:"
echo "${param_file}"
echo "Number of jobs: ${n_jobs}"
echo
