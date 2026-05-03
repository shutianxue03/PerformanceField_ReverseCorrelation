#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Submit OOD_sim Slurm array jobs in chunks.
#
# This avoids QOSMaxSubmitJobPerUserLimit when the full array is too large.
# ============================================================

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
param_file="${project_root}/Data_job_params/OOD_sim_params.tsv"
array_script="${project_root}/Codes/shell_Sim.sh"

# Settings
batch_size=500        # number of array tasks per submitted batch
max_running=50        # Slurm array throttle: %50
sleep_seconds=2       # small pause between sbatch calls

if [[ ! -f "${param_file}" ]]; then
  echo "Error: cannot find parameter file:"
  echo "  ${param_file}"
  echo "Run:"
  echo "  bash Codes/shell_Sim_make_param_table.sh"
  exit 1
fi

if [[ ! -f "${array_script}" ]]; then
  echo "Error: cannot find array script:"
  echo "  ${array_script}"
  exit 1
fi

# Count jobs: total lines minus header
n_jobs=$(( $(wc -l < "${param_file}") - 1 ))

if (( n_jobs <= 0 )); then
  echo "Error: parameter table has no jobs."
  exit 1
fi

echo "Parameter file: ${param_file}"
echo "Array script: ${array_script}"
echo "Total jobs: ${n_jobs}"
echo "Batch size: ${batch_size}"
echo "Max running per batch: ${max_running}"
echo "======================================"

start=1

while (( start <= n_jobs )); do
  end=$(( start + batch_size - 1 ))

  if (( end > n_jobs )); then
    end=${n_jobs}
  fi

  echo "Submitting array batch: ${start}-${end}%${max_running}"

  sbatch --array="${start}-${end}%${max_running}" "${array_script}"

  start=$(( end + 1 ))

  sleep "${sleep_seconds}"
done

echo "======================================"
echo "All array batches submitted."
echo "Total jobs requested: ${n_jobs}"
echo "======================================"