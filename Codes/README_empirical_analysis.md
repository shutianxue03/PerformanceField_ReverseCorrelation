# Empirical Data Analysis README (`shell_launch_NOM.sh`)

This document describes the trial-wise NOM analysis pipeline for **empirical (human) data**, launched by `shell_launch_NOM.sh`.

## What this pipeline does

The empirical NOM workflow runs in two stages:

1. **Stage 1 (`compIV`)**
  Compute trial-wise internal variables (DVs/IVs), split trials (template/train/test), derive templates, and save per-iteration data.
2. **Stage 2 (`fitNOM`)**
  Fit Model B parameters to stage-1 outputs and predict behavioral metrics.

Main launcher:

- `shell_launch_NOM.sh`

Wrappers called by launcher:

- `shell_NOM_compIV.sh`
- `shell_NOM_fitNOM.sh`

Core MATLAB functions:

- `OOD_NOM_Trialwise_compIV.m`
- `OOD_NOM_Trialwise_fitNOM.m`

---

## Before you run (important)

1. **Set environment paths** in `fxn_analysis_RC_v2/SX_RC1_setting.m`:
  - `nameFolder_server`
  - `str_envir` (`HPC`, `Server`, or `Local`)
2. **Verify empirical inputs already exist** for each subject:
  - behavioral file: `Data/Data_OOD_<nORI><nSF>/<subj><nblocks>/<subj>_behavMeas.mat`
  - energy file: `Data/Data_OOD_<nORI><nSF>/<subj><nblocks>/<subj>_energy_T_<nORI>_<nSF>.mat`
3. **Set stage switch** in `shell_launch_NOM.sh`:
  - `flag_stage=1` for compIV
  - `flag_stage=2` for fitNOM
4. **Interface check (current codebase):**
  - `OOD_NOM_Trialwise_compIV` expects  
   `(..., iLocComb, lambda_whiten, flag_regressType, iModelA, nIter, nJob, iJob)`
  - `OOD_NOM_Trialwise_fitNOM` expects  
  `(..., iLocComb, flag_incluCrit, C_contribution, iModelA, iModelB, nIter, nJob, iJob)`
  - Current `shell_NOM_compIV.sh` and `shell_NOM_fitNOM.sh` call older signatures. Update wrapper arguments before production runs.

---

## How to run on HPC (Slurm)

From `Codes/`:

1. Edit `shell_launch_NOM.sh` parameters:
  - `flag_stage`
  - `nIter`
  - `nJob`
  - `isubjList`
  - `iLocCombList`
  - `iModelAList`
  - `iModelBList` (used in stage 2)
2. Launch:

```bash
bash shell_launch_NOM.sh
```

The launcher submits one Slurm job per parameter combination with:

- stage 1: `sbatch shell_NOM_compIV.sh ...`
- stage 2: `sbatch shell_NOM_fitNOM.sh ...`

---

## Parameter conventions

### Subject index (`isubj`)

Mapped in MATLAB code (`OOD_NOM_Trialwise_compIV.m` / `fitNOM.m`) to:

`{'YK','SP','SX','LS','RE','MD','AS','HL','FH','HA','CS','DT','DU','RC','SR'}`

### Location index (`iLocComb`)

- `1..5`: single locations
- `6`: horizontal meridian (2+4)
- `7`: vertical meridian (3+5)
- `8`: perifovea (2:5)

(`shell_launch_NOM.sh` currently uses `1:5`.)

### Model A (`iModelA`)

From settings:

- `1`: RC-derived template
- `2`: ideal template
- `3`: randomized template (if enabled in specific scripts)

### Model B (`iModelB`)

In `fitNOM`, common variants include:

- `1`: full model
- `2`: no shared noise
- `3`: no multiplicative noise
- `4`: no additive noise
- `5/6/7`: reduced single-noise variants (available in MATLAB function)

---

## Output files

Outputs are written under:

- `Data/Data_NOM_Trialwise_<nORI><nSF>/`
- `Figures/Human/<subj>/`

Per subject/location:

- Stage 1 (`compIV`):
  - `n<nIter>_J<iJob>_A<iModelA>_compIV.mat`
  - contains split data, IVs, templates, criterion info, metrics across iterations
- Stage 2 (`fitNOM`):
  - `n<nIter>_J<iJob>_A<iModelA>B<iModelB>.mat`
  - contains fitted params, nLL, predicted metrics across iterations

Progress snapshots are also saved during loops with `_iIter_elapsedmin` suffixes.

---

## Typical execution order

1. Run **all stage-1 jobs** (`flag_stage=1`) and confirm `*_compIV.mat` exists.
2. Run **all stage-2 jobs** (`flag_stage=2`) using matching `nIter/nJob/iJob` settings.
3. Aggregate and plot with downstream analysis scripts (for example, scripts in `SX_analysis5_`* and plotting utilities under `fxn_NOM` / `fxn_analysis_RC_v2` depending on your figure target).

---

## Practical notes

- Use deterministic `nIter/nJob/iJob` grids consistently between stages.
- If running on HPC, plotting is mostly disabled in MATLAB functions by `str_envir`.
- Check Slurm logs (`zzz_NOM1_*.out`, `zzz_NOM2_*.out`) for missing-input errors first.
- Start with one subject/location/job as a smoke test before full grid submission.

# ===============================================================================

================================================================================

# Group Analysis README (`SX_analysis5_NOM_Trialwise.m`)

This README covers the stage where trial-wise NOM outputs are combined across observers, statistics are run, and group figures are generated.

## Main script

- `SX_analysis5_NOM_Trialwise.m`

This script is the post-fit analysis hub for empirical data. It assumes per-subject outputs from:

- `OOD_NOM_Trialwise_compIV.m`
- `OOD_NOM_Trialwise_fitNOM.m`

---

## What this script does

At a high level, it:

1. Defines subject cohort(s), model settings, location groups, and plotting/statistics options.
2. Loads a compiled multi-observer data file (`*_allCond` variables).
3. Appends observer-level behavioral summaries (`*_allSubj`) from `Data_OOD`.
4. Produces a full panel of figures (24 major sections in the script).
5. Runs repeated-measures model-comparison statistics on nLL.

---

## Critical prerequisite (read this first)

`SX_analysis5_NOM_Trialwise.m` **loads many `*_allCond` variables** from:

- `nameFolder_Data_SaveCompile = sprintf('%s/n%d_n%d', nameFolder_Data_NOM_Trialwise, nSubj, nIter*nJob)`

For example (with current defaults):

- `nIter=200`, `nJob=5` -> `nIter*nJob=1000`
- if `nSubj=12` -> expected file stem: `Data/Data_NOM_Trialwise_<nORI><nSF>/n12_n1000`

The large compile block that builds `*_allCond` in this script is currently commented out.  
So in normal use, you should run this script with an already prepared compile file.

If this file is missing, later `load(nameFolder_Data_SaveCompile, ...)` calls will fail.

---

## Inputs expected by this script

### 1) Compiled NOM arrays (`*_allCond`)

Examples loaded throughout the script:

- `metrics_allCond`
- `template_tmpl_allCond`, `template_full_allCond`
- `sep_allCond`
- `margORI_allCond`, `margSF_allCond`
- `margPred_*_allCond`, `margParams_*_allCond`, `margR2_*_allCond`
- `margTunC_*_allCond`
- `IV_allCond`, `nTrials_allCond`
- `metric_data_allCond`, `metric_pred_allCond`
- `R2_NOM_allCond`, `R2_w_NOM_allCond`
- `nLL_allCond`
- `params_allCond`

### 2) Behavioral source files (per observer)

Loaded from `Data_OOD` as:

- `.../<subj><nblocks>/<subj>_behavMeas.mat` (expects `*_perSess_perLoc` variables)

These are used to compute and append:

- `CS_allSubj`, `dprime_allSubj`, `criterion_allSubj`, `RT_allSubj`, `pC_allSubj`, `pA_allSubj`

---

## Configuration you will typically edit

Near the top of `SX_analysis5_NOM_Trialwise.m`:

- cohort/model:
  - `iRun` (selects subject inclusion set)
  - `iModelA_sim_all`, `iModelB_sim_all`
  - `iModelA_plot`, `iModelB_plot`, `iModelB_plot_all`
- resampling size:
  - `nIter`, `nJob`
- location sets:
  - `iLocGroups_all`, `iLocSingle_allSets`
- inferential settings:
  - `nPerm`, `nBoot`, `CI95`, `CI68`, seeds
- correlation mode:
  - `flag_UseRUseRho` (`useR` vs `useRho`)

Also verify environment/path in:

- `fxn_analysis_RC_v2/SX_RC1_setting.m`

---

## How to run

From `Codes/` in MATLAB:

```matlab
SX_analysis5_NOM_Trialwise
```

The script runs sequentially through plotting/stat sections (labeled `1/24` ... `24/24` in console prints).

---

## Stats included in this script

Primary model-comparison stats are in sections around `17/24` and `18/24`:

- Repeated-measures ANOVA on nLL via `rm2ANOVA_A1`
  - factors include Model B and Location (Model A fixed in those sections)
- Bootstrap confidence intervals for partial eta-squared (`eta2p`)
  - computed with `fxn_eta2p_rm2`
  - bootstrapped over observers

Additional correlation/permutation analyses are done through helper plotting/stat functions, e.g.:

- `basicFxn_drawBars_permutation`
- `basicFxn_drawCorr_permutation`
- `basicFxn_drawCorrAsym_permutation`
- `basicFxn_compAsym_permutation`

---

## Output locations

Base output folder:

- `Figures/NOM_Trialwise_n<nSubj>/`

Subfolders created by the script include:

- `Behav`
- `Template`
- `Separability`
- `TuningFxns_group`
- `TuningFxns_IDVD`
- `TuningCs`
- `NOMmetrics_A*`
- `nLL`
- `NOMparams`
- `Corr/...`
- `CorrAsym/...`
- `CompAsym/...`
- `pAdevCorr`

Compiled data file used/updated:

- `Data/Data_NOM_Trialwise_<nORI><nSF>/n<nSubj>_n<nIter*nJob>.mat`

---

## Practical tips

- First run with a small cohort (`iRun`) and reduced `nIter/nJob` to smoke-test.
- Keep `nIter`/`nJob` consistent with the files you compiled from stage-1/2 outputs.
- If you only need a subset of figures, temporarily guard sections with simple `if` flags to shorten runtime.
- Use `checkMissingFiles.m` before group analysis when troubleshooting incomplete model runs.

