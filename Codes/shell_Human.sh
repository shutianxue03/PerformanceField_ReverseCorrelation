#!/usr/bin/env bash
#SBATCH --job-name=OODhuman
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=Logs/Human_%A_%a.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=shutianxue30@gmail.com

# ============================================================
# Run one OOD_NOM_Human array task from the parameter table.
# ============================================================

set -euo pipefail

# -----------------------------
# Resolve project paths
# -----------------------------

if [[ -n "${SLURM_SUBMIT_DIR:-}" && -d "${SLURM_SUBMIT_DIR}/Codes" ]]; then
  project_root="${SLURM_SUBMIT_DIR}"
  script_dir="${project_root}/Codes"
else
  script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  project_root="$(cd "${script_dir}/.." && pwd)"
fi

param_file="${project_root}/Data_job_params/OOD_Human_params.tsv"
log_dir="${project_root}/Logs"

mkdir -p "${log_dir}"

# -----------------------------
# Safety checks
# -----------------------------

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "Error: SLURM_ARRAY_TASK_ID is not set." >&2
  exit 1
fi

if [[ ! -f "${param_file}" ]]; then
  echo "Error: cannot find parameter file:" >&2
  echo "  ${param_file}" >&2
  echo "Run the parameter-table script first, for example:" >&2
  echo "  bash Codes/shell_Human_make_param_table.sh" >&2
  exit 1
fi

n_jobs=$(( $(wc -l < "${param_file}") - 1 ))

if (( n_jobs <= 0 )); then
  echo "Error: parameter table has no jobs." >&2
  exit 1
fi

# +1 because line 1 is the header.
line_number=$((SLURM_ARRAY_TASK_ID + 1))
line="$(sed -n "${line_number}p" "${param_file}")"

if [[ -z "${line}" ]]; then
  echo "Error: no parameter row found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
  echo "Requested line number: ${line_number}" >&2
  echo "Parameter file: ${param_file}" >&2
  exit 1
fi

# -----------------------------
# Read parameters from table
# -----------------------------

IFS=$'\t' read -r isubj iLocComb flag_normDV nIter nJob iJob <<< "${line}"

job_name="S${isubj}_L${iLocComb}_norm${flag_normDV}_n${nIter}_J${iJob}of${nJob}"

# -----------------------------
# Print metadata
# -----------------------------

echo "======================================"
echo "Started: $(date)"
echo "Array job ID: ${SLURM_ARRAY_JOB_ID:-NA}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Job ID: ${SLURM_JOB_ID:-NA}"
echo "Job name: ${job_name}"
echo "Node: ${SLURMD_NODENAME:-NA}"
echo "Project root: ${project_root}"
echo "Script dir: ${script_dir}"
echo "Parameter file: ${param_file}"
echo "Total jobs: ${n_jobs}"
echo "Parameter row:"
echo "${line}"
echo "======================================"

# -----------------------------
# Run MATLAB
# -----------------------------

cd "${project_root}"

module load matlab/2025b

matlab -batch "try, addpath(genpath('Codes')); OOD_NOM_Human(${isubj}, ${iLocComb}, ${flag_normDV}, ${nIter}, ${nJob}, ${iJob}); catch ME, disp(getReport(ME,'extended')); exit(1); end; exit(0);"

echo "======================================"
echo "Finished: $(date)"
echo "Array task ${SLURM_ARRAY_TASK_ID} completed."
echo "======================================"