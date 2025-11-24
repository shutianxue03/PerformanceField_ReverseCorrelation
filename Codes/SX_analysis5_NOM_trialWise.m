% function SX_makeAllFigures_NOM()
% SX_makeAllFigures_NOM
% Master wrapper to generate all NOM figures.

close all; clc;
SX_RC1_setting;  %#ok<*NOPRT>

% ---- user settings ----
nIterations = 3;
iModelA_main = 1;        % core template
iModelA_IO   = 3;        % ideal observer template (for model comparison)
iModelB_full = 1;        % full dual-noise + rho model
iLocComb_all = [1 2];% HM, VM, LVM, UVM
subjList = {'YK','SP','SX','LS','RE','MD','AS','HL','FH','HA','CS','DT','DU','RC','SR'};
nSubj    = numel(subjList);

%% SECTION 1: Kernels & template recovery
SX_plot_NOM_kernels(nIterations, iModelA_main, iLocComb_all, subjList);

%% SECTION 2: Parameters vs location (ModelA=1, ModelB=1)
SX_plot_NOM_params_vs_loc(nIterations, iModelA_main, iModelB_full, iLocComb_all, subjList);

%% SECTION 3: Measured vs predicted metrics (pYES, pC, pA)
metricsToPlot = {'pYES','pC','pA'};
SX_plot_NOM_fit_metrics(nIterations, iModelA_main, iModelB_full, iLocComb_all, subjList, metricsToPlot);

%% SECTION 4: pA vs parameters
SX_plot_NOM_pA_vs_params(nIterations, iModelA_main, iModelB_full, iLocComb_all, subjList);

%% SECTION 5: Model comparison ΔBIC, ModelA=1 vs 3
iModelB_all = 1:7;
SX_plot_NOM_model_comparison(nIterations, iModelA_main, iModelA_IO, iModelB_all, iLocComb_all, subjList);

disp('All NOM figure templates executed.');

