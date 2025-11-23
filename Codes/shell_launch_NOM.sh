#!/bin/bash
# ============================================================
# OOD_NOM_Trialwise batch launcher
# Created by Shutian Xue on 11/23/2025
# Description: Submit NOM fitting jobs for all subjects and locations
# ============================================================

set -euo pipefail

iModelA=1          # 1 = core model (see OOD_NOM_Trialwise.m)
nIterations=20     # number of iterations per (subj, location, modelA)

for isubj in {1..15}; do
  for iLocComb in {1..8}; do
    echo "Submitting: subj=${isubj}, loc=${iLocComb}, modelA=${iModelA}, nIter=${nIterations}"
    sbatch shell_NOM.sh "${isubj}" "${iLocComb}" "${iModelA}" "${nIterations}"
  done
done