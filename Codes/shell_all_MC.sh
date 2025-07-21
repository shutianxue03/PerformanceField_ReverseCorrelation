#!/bin/bash
# Created by Shutian Xue on 07/21/2025
# Last modified by Shutian Xue on 07/21/2025

# OOD_MC_current(isubj, iiLoc_all_all, flag_PatchMode, flag_cutMapping, flag_mirrorMapping, ifeature, ifamily, MCmode)
flag_PatchMode=2
flag_cutMapping=0
flag_mirrorMapping=1
MCmode=1  # 1=10-fold CV, 2=LOOCV, 3=IC (AIC, AICc, BIC)

echo "Submitting jobs with: flag_PatchMode=$flag_PatchMode, flag_cutMapping=$flag_cutMapping, flag_mirrorMapping=$flag_mirrorMapping, MCmode=$MCmode"

for iSubj in $(seq 1 15); do
    for iiLoc_all_all in 1 2 3; do
        echo "Submitting: iSubj=$iSubj, iiLoc_all_all=$iiLoc_all_all"

        # ORI feature
        ifeature=1 # ORI
        for ifamily in 1 8; do # 1=Gaussian, 8=DoG
            echo "  ifeature=$ifeature (ORI), ifamily=$ifamily"
            sbatch shell_MC.sh $iSubj $iiLoc_all_all $flag_PatchMode $flag_cutMapping $flag_mirrorMapping $ifeature $ifamily $MCmode
        done

        # SF feature
        ifeature=2 # SF
        for ifamily in 2 3; do # 2=log parabola, 3=truncated log parabola
            echo "  ifeature=$ifeature (SF), ifamily=$ifamily"
            sbatch shell_MC.sh $iSubj $iiLoc_all_all $flag_PatchMode $flag_cutMapping $flag_mirrorMapping $ifeature $ifamily $MCmode
        done
    done
done