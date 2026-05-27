# Empirical Data Analysis README

This file summarizes the current human-data NOM pipeline.

## Main batch workflow

Start from the human parameter table:

- `shell_Human_make_param_table.sh`

Then submit the array jobs with:

- `shell_Human_submit_batches.sh`

Each array task runs:

- `shell_Human.sh`

which calls:

- `OOD_NOM_Human`

## What the pipeline does

For each subject and location condition, `OOD_NOM_Human` runs two stages:

1. `OOD_NOM_Trialwise_compDV_A12`
   - computes DVs
   - derives templates
   - saves trial-wise split data and template summaries
2. `OOD_NOM_Trialwise_fitNOM`
   - fits Model B variants
   - predicts behavioral metrics
   - saves fitted parameters and prediction outputs

## Human job table

The current parameter table is created by `shell_Human_make_param_table.sh` with columns:


## Subject and location conventions

### Subjects

`isubj` maps to:

`{'YK','SP','SX','LS','RE','MD','AS','HL','FH','HA','CS','DT','DU','RC','SR'}`

### Locations

- `1..5`: single locations
- `6`: horizontal meridian = `[2, 4]`
- `7`: vertical meridian = `[3, 5]`
- `8`: perifovea = `2:5`

## Before running

1. Check paths in `fxn_analysis_RC_v2/SX_RC1_setting.m`.
2. Make sure the human behavioral and energy files already exist.
3. Confirm the output folders are writable on local or HPC storage.

## Main outputs

For each subject and location, outputs are saved under the NOM trial-wise data folder.

- `n<nIter>_J<iJob>_A1_compDV.mat`
- `n<nIter>_J<iJob>_A2_compDV.mat`
- `n<nIter>_J<iJob>_A<iModelA>B<iModelB>.mat`

Figures are written under the human figure folder when plotting is enabled.

## Group analysis

Main script:

- `SX_analysis5_NOM_Trialwise.m`

This script loads the saved human NOM outputs across observers, compiles group-level arrays, and generates the group summary figures and statistics.

### Inputs expected by the script

The script expects compiled `*_allCond` arrays under the NOM trial-wise data folder. Common variables include:

- `metrics_allCond`
- `template_tmpl_allCond`, `template_full_allCond`
- `sep_allCond`
- `margORI_allCond`, `margSF_allCond`
- `margPred_*_allCond`, `margParams_*_allCond`, `margR2_*_allCond`
- `margTunC_*_allCond`
- `DV_allCond`, `nTrials_allCond`
- `metric_data_allCond`, `metric_pred_allCond`
- `R2_NOM_allCond`, `R2_w_NOM_allCond`
- `nLL_allCond`
- `params_allCond`

It also loads per-subject behavioral files from `Data_OOD` to build:

- `CS_allSubj`
- `dprime_allSubj`
- `criterion_allSubj`
- `RT_allSubj`
- `pC_allSubj`
- `pA_allSubj`

### Run it

From MATLAB:

```matlab
SX_analysis5_NOM_trialWise
```

### Practical notes

- Keep `nIter` and `nJob` consistent with the saved job outputs.
- The large compile block inside the script may still be commented out, so make sure the compiled data file already exists if you are not rebuilding it there.
- For troubleshooting incomplete runs, use `checkMissingFiles.m` first.

