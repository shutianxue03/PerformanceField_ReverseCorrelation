#!/bin/bash
# Created by Shutian Xue on 07/16/2025
# Last modified by Shutian Xue on 11/09/2025

# This script runs simulations for different models and parameters
for noiseCST in 0 0.1 0.2 0.3 0.4 0.5
do
    for gaborCST in 0.1 0.2 0.3 0.4 0.5
    do
        for nTrials in 2000 4000 8000
        do
        % write if statement on iModelB_sim to only run model 4
            if [ $iModelB_sim -eq 4 ]; then
                for noiseP in $(seq(0, 5, 20))
                do  
                    sbatch shell_Sim.sh $noiseCST $gaborCST $nTrials $noiseP $iModelB_sim
                done
            fi

            if [ $iModelB_sim -eq 5 ]; then
                for noiseP in $(seq(0, .1, .5))
                do  
                    sbatch shell_Sim.sh $noiseCST $gaborCST $nTrials $noiseP $iModelB_sim
                done
            fi
        done
    done    
done