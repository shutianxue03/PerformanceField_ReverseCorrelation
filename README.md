# PF_RC Project Overview

This repository contains the code and analysis pipeline for the PF_RC project: a combination of human psychophysical experimentation, model-based simulation, and NOM-based empirical data analysis. The project is organized around MATLAB scripts in `Codes/`, with supporting data folders under `Data_*` and figures under `Figures/`.

The repository is structured around three complementary workflows:

1. Human data collection via the psychophysical experiment pipeline.
2. Model simulation and recovery analysis.
3. Empirical NOM fitting and group-level analysis.

This document provides a single entry point for understanding how the pieces fit together.

---

## Project structure

- `Codes/` — core MATLAB scripts, batch shell runners, analysis functions, and project documentation.
- `Data_Human/` and `Data_Sim/` — human and simulated datasets.
- `Figures/` — plots and figure outputs.
- `Outputs/` — compiled result files and summary outputs.
- `Manuscript/` — manuscript and bibliography files.

The three README files in `Codes/` summarize the main components:

- [Codes/README_psychophysical_experiment.md](Codes/README_psychophysical_experiment.md)
- [Codes/README_model_simulation.md](Codes/README_model_simulation.md)
- [Codes/README_empirical_analysis.md](Codes/README_empirical_analysis.md)

---

## 1. Human psychophysical experiment pipeline

The experiment side is launched from:

- `PF_RC_runExp.m`

This script is the experiment entry point. It initializes the session, collects observer information, prepares experiment parameters, optionally initializes EyeLink, loads or creates the subject/session record, runs the blocks and trials, and then saves the session outcome and post-session summary.

### Typical run

From MATLAB in the `Codes/` directory:

```matlab
PF_RC_runExp
```

The program prompts for:

- subject ID / observer initials
- experiment mode
- EyeLink status
- demonstration mode flag
- threshold values when running the main experiment

### Main modes

- Titration: threshold estimation
- Experiment: main data collection
- Practice: training sessions
- Quick titration: short threshold update
- Simulation mode: special run path

### Recommended observer workflow

1. Run titration to estimate threshold.
2. Optionally run quick titration.
3. Run the main experiment using the selected threshold.
4. Use practice mode for training and demos as needed.

### Experiment details

- 5 possible stimulus locations
- Binary present/absent judgment task
- Response mapping on Mac: `F` = yes, `J` = no
- Titration and experiment block designs are parameterized by session and staircase settings

### Output files

Each subject has a data folder similar to:

- `Data/<subjName>/`

Files are named by run type and time stamp, such as:

- titration files
- experiment files
- practice files
- quick titration files

Additional outputs include:

- threshold posterior PDF files
- per-session parameter snapshots
- EyeLink data when enabled

This side of the project is intended to generate the behavioral measurements used downstream in the NOM analysis pipeline.

---

## 2. Simulation pipeline

The simulation side is centered on:

- `OOD_sim.m`

This script simulates the experiment under controlled conditions and then applies the same NOM-based downstream analysis pipeline used for the human data.

### Simulation entry point

```matlab
OOD_sim(noiseCST, gaborCST, nTrials, Nmul_true, Nadd_true, Nshared_true, cSDT_true, iModelB_sim, nIter)
```

### What the simulation does

`OOD_sim.m`:

1. simulates stimuli and orientation x spatial-frequency energy
2. builds the ground-truth template and decision variable
3. adds internal noise
4. generates behavioral responses and summary metrics
5. runs the same downstream NOM processing used for empirical data

### Example call

```matlab
OOD_sim(0.2, 0.3, 10000, 0.2, 1, 0.5, 0, 1, 50)
```

### Simulation parameters

- `noiseCST`: external noise contrast
- `gaborCST`: signal contrast
- `nTrials`: number of trials (must be even)
- `Nmul_true`, `Nadd_true`, `Nshared_true`: internal-noise parameters
- `cSDT_true`: criterion value in SDT units
- `iModelB_sim`: model family index
- `nIter`: resampling iterations

### Main simulation outputs

Outputs are organized into `Data_*` directories and include:

- truth files with the simulated template and parameters
- behavioral measurement files
- trial-wise NOM outputs saved as `*_A1_compDV.mat`, `*_A2_compDV.mat`, and model-fitting outputs

### Batch simulation workflow

Simulation batch parameters are generated using:

- `shell_Sim_make_param_table.sh`

and are submitted via the shell launch and batch scripts in `Codes/`.

These simulations are especially useful for testing model recovery, parameter sensitivity, and fit quality under known ground truth.

---

## 3. Empirical analysis and NOM fitting

The empirical analysis pipeline is organized around the human-data workflow.

### Main human job workflow

1. Create the parameter table:
   - `shell_Human_make_param_table.sh`
2. Submit jobs:
   - `shell_Human_submit_batches.sh`
3. Each task calls:
   - `shell_Human.sh`
4. The shell script runs:
   - `OOD_NOM_Human`

### What the analysis does

For each subject and location condition, `OOD_NOM_Human` runs:

1. `OOD_NOM_Trialwise_compDV_A12`
   - computes DVs
   - derives templates
   - saves trial-wise split data and template summaries
2. `OOD_NOM_Trialwise_fitNOM`
   - fits Model B variants
   - predicts behavioral metrics
   - saves fitted parameters and prediction outputs

### Subject and location conventions

The project tracks human subjects using the `isubj` mapping, including observers such as YK, SP, SX, LS, RE, MD, AS, HL, FH, HA, CS, DT, DU, RC, and SR.

Locations are encoded as:

- `1..5`: single locations
- `6`: horizontal meridian
- `7`: vertical meridian
- `8`: perifovea

### Group analysis

The primary group-level analysis script is:

- `SX_analysis5_NOM_Trialwise.m`

This script loads outputs across subjects and conditions, compiles arrays, and generates group summary figures and statistics.

It expects compiled arrays such as:

- `metrics_allCond`
- `template_tmpl_allCond`
- `template_full_allCond`
- `sep_allCond`
- `margORI_allCond`
- `margSF_allCond`
- `DV_allCond`
- `params_allCond`
- `R2_NOM_allCond`
- `nLL_allCond`

and also combines per-subject behavioral summaries such as:

- `CS_allSubj`
- `dprime_allSubj`
- `criterion_allSubj`
- `RT_allSubj`
- `pC_allSubj`
- `pA_allSubj`

### Run group analysis

From MATLAB:

```matlab
SX_analysis5_NOM_trialWise
```

---

## Typical project workflow

A typical full project cycle is:

1. Collect human psychophysical data with `PF_RC_runExp.m`.
2. Run or prepare simulation conditions with `OOD_sim.m`.
3. Fit NOM models on human and/or simulated data.
4. Compile outputs across subjects and conditions.
5. Run group analyses and generate summary figures.

This provides both a behavioral study and a model-recovery framework within the same repository.

---

## Setup and practical notes

Before running the pipeline:

- confirm MATLAB paths and project settings in `fxn_analysis_RC_v2/SX_RC1_setting.m`
- ensure the relevant `Data_*` folders exist and are writable
- verify that output folders are available on local or HPC storage
- check that the right number of iterations and job IDs match the saved files

### Useful maintenance scripts

- `checkMissingFiles.m` — useful for troubleshooting incomplete runs
- plotting and summary scripts under `Codes/` when generating figures and diagnostics

### Notes on output conventions

The project stores internally generated files under condition-specific IO folders and uses a consistent trial-wise naming convention across simulation and empirical analyses. The downstream fitting stage operates on saved DV outputs rather than older IV-based files.

---

## Summary

This project bridges three tightly related workflows:

- collecting human psychophysical observer data,
- simulating the same task under controlled parameter settings,
- and fitting NOM-based models to summarize behavior and compare model explanations across conditions.

Together, these workflows support a full analysis pipeline for perceptual decision-making, internal noise, template-based coding, and model recovery.
