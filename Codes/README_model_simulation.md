# Model Simulation README (`OOD_sim.m`)

This document covers the model-simulation pipeline driven by `OOD_sim.m`.

## What this pipeline does

`OOD_sim.m` simulates trial-wise observer behavior using a noisy observer model (NOM), then runs the same analysis/fitting pipeline used for empirical data:

1. Generate simulated stimuli and ORI×SF energy.
2. Build clean decision variables (DV) from a ground-truth template.
3. Add internal noise (multiplicative/additive/shared).
4. Set criterion to match a target criterion in z-space (`Cz_true`).
5. Generate binary responses and behavioral metrics.
6. Run:
   - `OOD_NOM_Trialwise_compIV` (template + IV computation),
   - `OOD_NOM_Trialwise_fitNOM` (parameter fitting + predicted metrics).

---

## Main entry point

- `OOD_sim.m`

```matlab
OOD_sim(noiseCST, gaborCST, nTrials, ...
    Nmul_true, Nadd_true, Nshared_true, Cz_true, ...
    lambda_whiten, flag_regressType, flag_incluCrit, C_contribution, ...
    iModelA_sim, iModelB_sim, nIter)
```

---

## Input arguments (recommended interpretation)

- `noiseCST`: external noise contrast (0-1).
- `gaborCST`: signal (Gabor) contrast (0-1).
- `nTrials`: total number of trials (must be even; code uses pair structure).
- `Nmul_true`: multiplicative internal noise scale.
- `Nadd_true`: additive internal noise SD.
- `Nshared_true`: shared internal noise SD across passes.
- `Cz_true`: criterion in SDT z-units (target criterion).
- `lambda_whiten`: whitening shrinkage used in `OOD_NOM_Trialwise_compIV`.
- `flag_regressType`: regression mode in compIV.
  - `1`: univariate
  - `2`: multivariate + smoothing
- `flag_incluCrit`: whether criterion is included as a fitted NOM parameter.
  - `0`: fit noise params only
  - `1`: fit noise params + criterion
- `C_contribution`: weight for criterion mismatch term in objective.
- `iModelA_sim`: template mode for fitting stage.
  - `1`: RC-derived template
  - `2`: ideal (true) template
  - (`3` exists in settings as randomized template)
- `iModelB_sim`: NOM variant for both simulation/fitting.
  - `1`: full (multi + add + shared)
  - `2`: no shared
  - `3`: no multiplicative
  - `4`: no additive
  - (fit stage also supports reduced models 5/6/7)
- `nIter`: number of resampling iterations in compIV/fitNOM.

---

## Required setup before running

1. Open MATLAB in the `Codes` directory.
2. Check environment paths in `fxn_analysis_RC_v2/SX_RC1_setting.m`.
   - This file defines `nameFolder_server` and `str_envir`.
   - For local runs, make sure local path is active.
3. Ensure target output directories are writable (they are auto-created if missing):
   - `Data/Data_OOD_<nORI><nSF>/...`
   - `Data/Data_NOM_Trialwise_<nORI><nSF>/...`
   - `Figures/IO/...`
4. If using parallel loops, start a MATLAB parallel pool if desired.

---

## Minimal run example

From `Codes/`:

```matlab
OOD_sim( ...
    0.2, ...      % noiseCST
    0.2, ...      % gaborCST
    4000, ...     % nTrials
    0.1, ...      % Nmul_true
    10, ...       % Nadd_true
    0, ...        % Nshared_true
    -0.5, ...     % Cz_true
    0.5, ...      % lambda_whiten
    1, ...        % flag_regressType
    1, ...        % flag_incluCrit
    1.0, ...      % C_contribution
    1, ...        % iModelA_sim
    1, ...        % iModelB_sim
    20);          % nIter
```

---

## Key outputs

For each simulation condition, `OOD_sim.m` creates an IO name like:

`IO_cNxx_cGxx_nT..._Nm..._Na..._Ns..._Cz..._cont..._whiten..._R..._B...`

and writes to:

- `Data_OOD.../<nameIO>/`
  - `truth.mat` (true template + criterion in DV units)
  - `behavMeas.mat` (`dataMatrix`, `metrics_sim`)
  - temporary `energy_T_<nORI>_<nSF>.mat` (deleted at end to save space)
- `Data_NOM_Trialwise.../<nameIO>/`
  - `n<nIter>_J<iJob>_A<iModelA>_compIV.mat`
  - `n<nIter>_J<iJob>_A<iModelA>B<iModelB>.mat` (fit results)
- `Figures/IO/<nameIO>/`
  - distribution/performance and fit summary figures (unless suppressed on HPC)

---

## Data conventions used in simulation

- Two-pass paired design:
  - `nPairs = nTrials / 2`
  - pass A is simulated, then copied to pass B at stimulus/energy level
  - shared noise creates across-pass correlation in noisy DV
- `dataMatrix` follows project conventions (trial, session, PRS/ABS, pass, pair, response, contrast, etc.).
- Behavioral metrics saved in `metrics_sim`:
  - `[dprime, criterion_z, pC, pHit, pFA, pA, pYES]`

---

## Pipeline dependencies

- Simulation + stimulus:
  - `exp_CreateFilteredNoise`, `exp_CreateGabor`, `exp_CreateCircularApertureSin`
- Energy/filtering:
  - `SX_sim02_setFilters`, `SX_RC4_Energy_parfor`
- IV/template:
  - `fxn_getIV_v3`, `fxn_getTemplate`
- SDT/metrics:
  - `SX_sim06_SDT`
- Downstream fitting:
  - `OOD_NOM_Trialwise_compIV`
  - `OOD_NOM_Trialwise_fitNOM`

---

## Notes and practical tips

- `nTrials` must be even.
- Keep `rng(1)` for reproducible simulation unless you intentionally vary seeds.
- If running on HPC, plotting is automatically suppressed in several scripts.
- Model B consistency: simulation currently fits with the same `iModelB_sim` by default.
- Temporary energy files are deleted at the end of `OOD_sim.m`; if debugging compIV inputs, comment out that deletion step.

