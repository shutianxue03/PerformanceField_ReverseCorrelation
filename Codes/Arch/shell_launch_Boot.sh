#!/bin/bash
# Created by Shutian Xue on 07/21/2025
# Last modified by Shutian Xue on 07/21/2025

# OOD_boot_current(iBatch, iSubj, iiLoc_all_all, nB_perBatch, ifamilyORI, ifamilySF, flag_cutMapping, flag_mirrorMapping)
nB_perBatch=100
ifamilyORI=1
ifamilySF=2
flag_cutMapping=0
flag_mirrorMapping=1

for iBatch in $(seq 1 10); do
    for iSubj in $(seq 1 15); do
        for iiLoc_all_all in 1 2 3; do
            echo "Submitting: iBatch=$iBatch, iSubj=$iSubj, iiLoc_all_all=$iiLoc_all_all"
            sbatch shell_Boot.sh $iBatch $iSubj $iiLoc_all_all $nB_perBatch $ifamilyORI $ifamilySF $flag_cutMapping $flag_mirrorMapping
        done
    done
done