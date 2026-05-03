#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Check completion status of all OOD_sim Slurm array jobs.
#
# This script automatically detects all Slurm array jobs matching
# a job-name pattern, summarizes their statuses, and writes one
# combined report.
#
# Default behavior:
#   - job name pattern: OODsim
#   - start time: today
#
# Usage:
#   bash Codes/shell_check_Sim_all_arrays.sh
#
# Optional:
#   bash Codes/shell_check_Sim_all_arrays.sh "2026-05-02"
#   bash Codes/shell_check_Sim_all_arrays.sh "2026-05-02T00:00:00"
#
# If your Slurm job name is not OODsim, edit job_name_pattern below.
# ============================================================

# -----------------------------
# User settings
# -----------------------------

job_name_pattern="OODsim"

# If provided, use user-specified start time.
# Otherwise, use today.
start_time="${1:-today}"

# -----------------------------
# Resolve paths
# -----------------------------

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "${script_dir}/.." && pwd)"

param_file="${project_root}/Data_job_params/OOD_sim_params.tsv"
summary_dir="${project_root}/Logs"
summary_file="${summary_dir}/OODsim_all_arrays_summary_$(date +%Y%m%d_%H%M%S).txt"

mkdir -p "${summary_dir}"

# -----------------------------
# Expected job count from table
# -----------------------------

if [[ -f "${param_file}" ]]; then
  n_expected=$(( $(wc -l < "${param_file}") - 1 ))
else
  n_expected="NA"
fi

# -----------------------------
# Temporary files
# -----------------------------

tmp_sacct="$(mktemp)"
tmp_squeue="$(mktemp)"
tmp_combined="$(mktemp)"
tmp_parent_ids="$(mktemp)"

cleanup() {
  rm -f "${tmp_sacct}" "${tmp_squeue}" "${tmp_combined}" "${tmp_parent_ids}"
}
trap cleanup EXIT

# -----------------------------
# Query sacct
# -----------------------------
# -X removes batch/extern substeps.
# --parsable2 gives pipe-delimited output.
#
# We query by user + job name + start time.
# This should capture all completed/failed/running array tasks
# that Slurm accounting knows about.

sacct \
  -u "$USER" \
  --name="${job_name_pattern}" \
  --starttime="${start_time}" \
  -X \
  --parsable2 \
  --noheader \
  --format=JobIDRaw,JobID,JobName%40,State,ExitCode,Elapsed,MaxRSS,Submit,Start,End \
  > "${tmp_sacct}" || true

# -----------------------------
# Query squeue
# -----------------------------
# -r expands array jobs so each array task gets one row.
# This is important because otherwise Slurm may show compressed ranges.
#
# Format:
# parentID|jobID|jobName|state|elapsed|reason

squeue \
  -u "$USER" \
  -r \
  -h \
  -o "%A|%i|%j|%T|%M|%R" \
  > "${tmp_squeue}" || true

# Keep only matching job names.
awk -F'|' -v pat="${job_name_pattern}" '$3 ~ pat {print $0}' "${tmp_squeue}" > "${tmp_squeue}.filtered"
mv "${tmp_squeue}.filtered" "${tmp_squeue}"

# -----------------------------
# Detect parent array IDs
# -----------------------------

# From sacct array-task IDs, e.g. 7823070_1 -> 7823070
awk -F'|' '
  $2 ~ /_[0-9]+$/ {
    split($2, a, "_");
    print a[1];
  }
' "${tmp_sacct}" >> "${tmp_parent_ids}"

# From squeue parent ID column.
awk -F'|' '{print $1}' "${tmp_squeue}" >> "${tmp_parent_ids}"

sort -u "${tmp_parent_ids}" -o "${tmp_parent_ids}"

n_parent_arrays=$(wc -l < "${tmp_parent_ids}" | tr -d ' ')

# -----------------------------
# Build combined task-status table
# -----------------------------
# Columns:
# source | parentID | taskID | jobID | jobName | state | exitCode | elapsed | maxRSS | extra
#
# Priority:
#   squeue active states should override sacct, because squeue is live.
#   Then sacct fills in completed/failed/etc.
#
# We create both, then use awk to keep one row per jobID,
# preferring squeue rows over sacct rows.

# sacct rows
awk -F'|' '
  $2 ~ /_[0-9]+$/ {
    split($2, a, "_");
    parent=a[1];
    task=a[2];
    source="sacct";
    jobid=$2;
    jobname=$3;
    state=$4;
    exitcode=$5;
    elapsed=$6;
    maxrss=$7;
    extra="submit=" $8 ",start=" $9 ",end=" $10;
    print source "|" parent "|" task "|" jobid "|" jobname "|" state "|" exitcode "|" elapsed "|" maxrss "|" extra;
  }
' "${tmp_sacct}" > "${tmp_combined}.sacct"

# squeue rows
awk -F'|' '
  {
    parent=$1;
    jobid=$2;
    jobname=$3;
    state=$4;
    elapsed=$5;
    reason=$6;

    task="NA";
    if (jobid ~ /_[0-9]+$/) {
      split(jobid, a, "_");
      task=a[2];
    }

    source="squeue";
    exitcode="NA";
    maxrss="NA";
    extra="reason=" reason;

    print source "|" parent "|" task "|" jobid "|" jobname "|" state "|" exitcode "|" elapsed "|" maxrss "|" extra;
  }
' "${tmp_squeue}" > "${tmp_combined}.squeue"

cat "${tmp_combined}.squeue" "${tmp_combined}.sacct" | \
awk -F'|' '
  # Because squeue rows are listed first, keep first occurrence per jobID.
  !seen[$4] {
    seen[$4]=1;
    print $0;
  }
' > "${tmp_combined}"

rm -f "${tmp_combined}.sacct" "${tmp_combined}.squeue"

# -----------------------------
# Count statuses
# -----------------------------

n_detected=$(wc -l < "${tmp_combined}" | tr -d ' ')

n_completed=$(awk -F'|' '$6 ~ /^COMPLETED/ {n++} END {print n+0}' "${tmp_combined}")
n_running=$(awk -F'|' '$6 ~ /^RUNNING/ {n++} END {print n+0}' "${tmp_combined}")
n_pending=$(awk -F'|' '$6 ~ /^PENDING/ {n++} END {print n+0}' "${tmp_combined}")
n_failed=$(awk -F'|' '$6 ~ /^FAILED/ {n++} END {print n+0}' "${tmp_combined}")
n_cancelled=$(awk -F'|' '$6 ~ /^CANCELLED/ {n++} END {print n+0}' "${tmp_combined}")
n_timeout=$(awk -F'|' '$6 ~ /^TIMEOUT/ {n++} END {print n+0}' "${tmp_combined}")
n_oom=$(awk -F'|' '$6 ~ /^OUT_OF_MEMORY/ {n++} END {print n+0}' "${tmp_combined}")
n_nodefail=$(awk -F'|' '$6 ~ /^NODE_FAIL/ {n++} END {print n+0}' "${tmp_combined}")

n_active=$(( n_running + n_pending ))
n_unsuccessful=$(( n_detected - n_completed - n_running - n_pending ))

if [[ "${n_expected}" != "NA" ]]; then
  n_missing=$(( n_expected - n_detected ))
else
  n_missing="NA"
fi

# -----------------------------
# Write summary
# -----------------------------

{
  echo "======================================"
  echo "OOD_sim all-array completion summary"
  echo "Generated: $(date)"
  echo "User: ${USER}"
  echo "Project root: ${project_root}"
  echo "Job name pattern: ${job_name_pattern}"
  echo "Start time filter: ${start_time}"
  echo "Parameter file: ${param_file}"
  echo "======================================"
  echo

  echo "Detected parent array job IDs:"
  if (( n_parent_arrays > 0 )); then
    cat "${tmp_parent_ids}"
  else
    echo "None detected"
  fi
  echo

  echo "Counts:"
  echo "  Expected jobs from parameter table = ${n_expected}"
  echo "  Detected array tasks              = ${n_detected}"
  echo "  Missing from Slurm query          = ${n_missing}"
  echo
  echo "Status counts:"
  echo "  COMPLETED                         = ${n_completed}"
  echo "  RUNNING                           = ${n_running}"
  echo "  PENDING                           = ${n_pending}"
  echo "  FAILED                            = ${n_failed}"
  echo "  CANCELLED                         = ${n_cancelled}"
  echo "  TIMEOUT                           = ${n_timeout}"
  echo "  OUT_OF_MEMORY                     = ${n_oom}"
  echo "  NODE_FAIL                         = ${n_nodefail}"
  echo "  UNSUCCESSFUL finished/known jobs  = ${n_unsuccessful}"
  echo

  if [[ "${n_expected}" != "NA" ]]; then
    if (( n_completed == n_expected )); then
      echo "Result: ALL EXPECTED JOBS COMPLETED SUCCESSFULLY."
    elif (( n_active > 0 )); then
      echo "Result: JOBS ARE STILL RUNNING OR PENDING."
    elif (( n_detected < n_expected )); then
      echo "Result: SOME EXPECTED JOBS WERE NOT DETECTED."
      echo "Interpretation: they may not have been submitted, may be outside the start-time filter, or sacct has not indexed them yet."
    else
      echo "Result: ALL DETECTED JOBS ARE FINISHED, BUT NOT ALL COMPLETED SUCCESSFULLY."
    fi
  else
    if (( n_active > 0 )); then
      echo "Result: JOBS ARE STILL RUNNING OR PENDING."
    else
      echo "Result: NO ACTIVE JOBS FOUND."
    fi
  fi

  echo
  echo "======================================"
  echo "Non-completed / active task list"
  echo "source|parentID|taskID|jobID|jobName|state|exitCode|elapsed|maxRSS|extra"
  echo "======================================"
  awk -F'|' '$6 !~ /^COMPLETED/ {print $0}' "${tmp_combined}"

  echo
  echo "======================================"
  echo "Failed/cancelled/timeout/OOM task list"
  echo "source|parentID|taskID|jobID|jobName|state|exitCode|elapsed|maxRSS|extra"
  echo "======================================"
  awk -F'|' '$6 ~ /^(FAILED|CANCELLED|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL)/ {print $0}' "${tmp_combined}"

  echo
  echo "======================================"
  echo "All detected task records"
  echo "source|parentID|taskID|jobID|jobName|state|exitCode|elapsed|maxRSS|extra"
  echo "======================================"
  cat "${tmp_combined}"

} | tee "${summary_file}"

echo
echo "Summary saved to:"
echo "${summary_file}"