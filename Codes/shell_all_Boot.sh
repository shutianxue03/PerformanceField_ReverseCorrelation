#!/bin/bash
# Created by Shutian Xue on 07/21/2025
# Last modified by Shutian Xue on 07/21/2025

# OOD_boot_current(ibatch, isubj, iiLoc_all_all, nB_perBatch, ifamilyORI, ifamilySF, flag_cutMapping, flag_mirrorMapping)
nB_perBatch=100
ifamilyORI=1
ifamilySF=2
flag_cutMapping=0
flag_mirrorMapping=1

for ibatch in 0 0.1 0.2 0.3 0.4 0.5; do
    for isubj in 0.1 0.2 0.3 0.4 0.5; do
        for iiLoc_all_all in 1000 5000 10000; do
            sbatch shell_Boot.sh $ibatch $isubj $iiLoc_all_all $nB_perBatch $ifamilyORI $ifamilySF $flag_cutMapping $flag_mirrorMapping
        done
    done
done