# Psychophysical Experiment README (`PF_RC_runExp.m`)

This README is for running the **human psychophysical experiment** pipeline using:

- `PF_RC_runExp.m`

It summarizes setup, run modes, expected prompts, output files, and practical checks.

---

## What `PF_RC_runExp.m` does

`PF_RC_runExp.m` is the experiment entry point. It:

1. Initializes paths and Psychtoolbox.
2. Collects operator inputs (subject ID, mode, eye-tracker flag).
3. Builds session parameters (`SX1_initParams`).
4. Optionally initializes EyeLink (`SX3_initEL`).
5. Loads/creates session record (`SX4_checkFile`).
6. Runs all blocks/trials (`SX7_run_AllBlocks`).
7. Ends session and closes devices (`SX12_endExp`).
8. Performs post-session save/cleanup and titration plotting (`SX13_postEXP`).

---

## Requirements

- MATLAB with Psychtoolbox (project notes mention PTB `3.0.18`).
- Local repository structure intact under `Codes/`.
- Working `Data/` directory (script auto-creates subject folder if needed).
- If eye tracking is enabled, EyeLink setup and write access to `eyedata/`.

---

## Run instructions

From MATLAB in `Codes/`:

```matlab
PF_RC_runExp
```

You will be prompted for:

- `subjName` (observer initials)
- `exp_mode`
  - `0`: Titration
  - `1`: Experiment
  - `2`: Practice
  - `3`: Quick titration
  - `4`: Simulation (listed in prompt; behavior differs from human run path)
- `EL_mode`
  - `1`: Eye-tracker ON
  - `0`: Eye-tracker OFF
- `demoFlag` (practice only)

For `exp_mode = 1` (main experiment), it also prompts:

- `thresh_exp` (threshold value to use in experiment blocks)

---

## Recommended session order (human observer)

1. Run **Titration** (`exp_mode=0`) for threshold estimation.
2. Optionally run **Quick titration** (`exp_mode=3`) at session start.
3. Run **Experiment** (`exp_mode=1`) with selected threshold.
4. Use **Practice** (`exp_mode=2`) for training/demo as needed.

---

## Experiment design highlights (from current code)

- 5 possible locations (fovea + 4 perifoveal).
- Task: binary present/absent (yes/no).
- Key mapping (Mac): `F` = yes, `J` = no.
- Main defaults:
  - Titration: 50 trials/stair/location, 2 staircases
  - Experiment: 100 trials/stair/location, 1 staircase
  - Practice: shorter blocks, no persistent data save
  - Quick titration: short staircase update
- `SX7_run_AllBlocks` requeues invalid fixation trials to block end.

---

## Output files and naming

Base folder per observer:

- `Data/<subjName>/`

Per-block/session files are named by mode prefix:

- Titration: `<subj>_stair_B###L#_yyyymmddTHHMM.mat`
- Experiment: `<subj>_exp_B###L#_yyyymmddTHHMM.mat`
- Practice: `<subj>_practice_B###L#_yyyymmddTHHMM.mat` (runtime mode exists; practice save is minimized)
- Quick titration: `<subj>_qstair_B###L#_yyyymmddTHHMM.mat`

Additional files:

- Threshold posterior PDF:
  - `Data/<subj>/<subj>_PDF.mat`
  - generated/updated by titration plotting utilities
- Session params snapshot:
  - `Data/<subj>/<subj>params<exp_mode>.mat`
- Eye data (if enabled):
  - renamed to match run file stem under `eyedata/`

---

## Important mode-specific behavior

- **Titration (`0`)**
  - If previous titration exists, it resumes and can reuse saved posterior.
  - At post-step, redundant intermediate titration files can be consolidated.
- **Experiment (`1`)**
  - Uses provided `thresh_exp`.
  - Creates stimuli in advance per block (`exp_createStim4AllTrials`).
- **Practice (`2`)**
  - Intended for training; data persistence is limited.
- **Quick titration (`3`)**
  - Loads prior PDF (`<subj>_PDF`) and updates quickly.

---

## Troubleshooting / sanity checks

- If `Data/<subj>/` does not exist, first titration run should create it.
- If quick titration fails loading PDF, run full titration first.
- If eye tracker is enabled and fails, rerun with `EL_mode=0` to isolate experiment logic.
- The script sets `Screen('Preference','SkipSyncTests',2)`; use stricter sync settings for final timing-critical collections if your lab machine is stable.
- Confirm key mappings on your OS (`SX1_initParams` handles Mac/Linux key codes).

---

## Related scripts (experiment side)

- `PF_RC_runExp.m` (entry point)
- `fxn_exp/SX1_initParams.m` (design/timing/staircases/keys)
- `fxn_exp/SX4_checkFile.m` (resume/create records)
- `fxn_exp/SX7_run_AllBlocks.m` (block/trial loop)
- `fxn_exp/SX11_endBlock.m` (block summary + save)
- `fxn_exp/SX12_endExp.m` (device shutdown + summary print)
- `fxn_exp/SX13_postEXP.m` (post-save, PDF plotting, eye-file rename)

