#!/bin/bash
# Created by Shutian Xue on 07/16/2025
# Last modified by Shutian Xue on 07/16/2025

for noiseCST in 0 0.1 0.2 0.3 0.4 0.5
do
    for gaborCST in 0.1 0.2 0.3 0.4 0.5
    do
        for nTrials in 1000 5000 10000
        do
            for noiseP in 0.1 0.2 0.3 0.4 0.5
            do  
                for iModelB_sim in 4
                do
                    sbatch shell_Sim.sh $noiseCST $gaborCST $nTrials $noiseP $iModelB_sim
                done
            done
        done
    done    
done