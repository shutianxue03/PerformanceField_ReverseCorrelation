#!/bin/bash
# ============================================================
# Created by Shutian Xue on 07/16/2025
# Last modified on 11/09/2025
# Description: Run simulations for different models and parameters
# ============================================================

for noiseCST in 0.2; do
  for gaborCST in $(seq 0.2 0.1 0.4); do
    for nTrials in 8000; do
      
      # ---- Run only model 4 ----
      ModelB_sim=4
      for noiseP in $(seq 0 5 20); do
        sbatch shell_Sim.sh "$noiseCST" "$gaborCST" "$nTrials" "$noiseP" "$ModelB_sim"
      done

      # ---- Run model 5 ----
      ModelB_sim=5
      for noiseP in $(seq 0 0.1 0.4); do
        sbatch shell_Sim.sh "$noiseCST" "$gaborCST" "$nTrials" "$noiseP" "$ModelB_sim"
      done

    done
  done
done
