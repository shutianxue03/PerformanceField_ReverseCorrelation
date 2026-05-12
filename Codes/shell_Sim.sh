#!/usr/bin/env bash
#SBATCH --job-name=OODsim
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=24G
#SBATCH --time=3:30:00
#SBATCH --output=Logs/Sim_%A_%a.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=shutianxue30@gmail.com

# ============================================================
# Created by Shutian Xue on 07/16/2025
# Last modified on 05/02/2026
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

param_file="${project_root}/Data_job_params/OOD_sim_params.tsv"
log_dir="${project_root}/Logs"

mkdir -p "${log_dir}"

# -----------------------------
# Safety checks
# -----------------------------

if [[ ! -f "${param_file}" ]]; then
  echo "Error: cannot find parameter file:" >&2
  echo "  ${param_file}" >&2
  echo "Run the parameter-table script first, for example:" >&2
  echo "  bash Codes/shell_Sim_make_param_table.sh" >&2
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

IFS=$'\t' read -r \
  noiseCST \
  gaborCST \
  nTrials \
  Nmul_true \
  Nadd_true \
  Nshared_true \
  cSDT_true \
  lambda_whiten \
  flag_regressType \
  flag_incluCrit \
  C_contribution \
  iModelB_sim \
  nIter \
  nBasisORI \
  nBasisSF <<< "${line}"

job_name="gCST${gaborCST}_Nm${Nmul_true}_Na${Nadd_true}_Ns${Nshared_true}_Cz${cSDT_true}_Bsim${iModelB_sim}"

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
echo "Parameter row:"
echo "${line}"
echo "======================================"

# -----------------------------
# Run MATLAB
# -----------------------------

cd "${script_dir}"

module load matlab/2025b

matlab -batch "try, OOD_sim(${noiseCST}, ${gaborCST}, ${nTrials}, ${Nmul_true}, ${Nadd_true}, ${Nshared_true}, ${cSDT_true}, ${lambda_whiten}, ${flag_regressType}, ${flag_incluCrit}, ${C_contribution}, ${iModelB_sim}, ${nIter}, ${nBasisORI}, ${nBasisSF}); catch ME, disp(getReport(ME,'extended')); exit(1); end; exit(0);"

echo "======================================"
echo "Finished: $(date)"
echo "Array task ${SLURM_ARRAY_TASK_ID} completed."
echo "======================================"
