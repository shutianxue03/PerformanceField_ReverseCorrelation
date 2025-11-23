function OOD_NOM_Trialwise(isubj, iLocComb, iModelA, nIterations)
% OOD_NOM_trialWise: This function generates IVs and fits different models to them.
% Inputs:
%      isubj: index of subj
%      iLocComb: 1=Fovea, 8=periF(6 deg ecc), 6=HM, 7=VM, 5=LVM, 3=UVM
%      iModelA: 1=core model, 2=permuted template, 3=use IO template
%      ni: number of iterations

clc
rng(123);   % define see for reproducibility
time_start = datetime('now');

%% Generate IVs
fprintf('\n\n=====================\nGenerating IVs (Model A%d)\n=====================\n\n', iModelA)

% The settings below can be placed in SX_RC1_setting
% IVType = 1;            % 1=sum of the dot product/convolution; 2=max; 3=normalized
% templateType = 3; % (1) raw (2) reconstructed kernel (3) mirrored template
% itype_template = 2; % 1=estimate template from PRS trials, ABS trials, or BOTH trials
% flag_PatchMode = 1; % if flag_PatchMode == 1, patchMode = 'T'; else, patchMode = 'N'; end

%------------------------------%
OOD_NOM_Trialwise_compIV(isubj, iLocComb, iModelA, nIterations);
%------------------------------%

%% Fit different models to IVs and make predictions
iModelB_all = 1:7; % Check the full list of models in SX_RC1_setting, search "namesModelB"
% iModelB_all = 2;
for iModelB = iModelB_all
    fprintf('\n\n=====================\nModel B%d\n=====================\n', iModelB)
    % ------------------------------%
    OOD_NOM_Trialwise_fitNOM(isubj, iLocComb, iModelA, iModelB, nIterations)
    % ------------------------------%
end

%% Ending
time_end = datetime('now');
elapsed = time_end - time_start;
fprintf('\n\n*** TOTAL time used: %s)\n', char(elapsed));
