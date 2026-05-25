function OOD_NOM_Human(isubj, iLocComb, lambda_whiten, nIter, nJob, iJob) %#ok<INUSD>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OOD_NOM_Human.m
%
% Human-data driver for the trial-wise NOM pipeline.
% 1) Compute DVs/templates via OOD_NOM_Trialwise_compDV_A12 (A1+A2)
% 2) Fit Model B variants for each Model A (A=1:2, B=1:7)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all;
warning off;
format compact;

time_start = datetime('now');
fprintf('\n\n%s: OOD_NOM_Human starts\n\n', datetime('now'))

% Add paths for custom functions.
addpath(genpath('Codes/'));

% Global settings (names/folders/flags shared across pipeline stages).
SX_RC1_setting;

% Defensive fallbacks for environments where script-defined vars are not visible to static checks.
if ~exist('namesLocComb', 'var'), namesLocComb = {}; end
if ~exist('C_contribution', 'var'), C_contribution = nan; end
if ~exist('flag_regressType', 'var'), flag_regressType = nan; end

% Model sets to run.
iModelA_fit_all = 1:2;
iModelB_fit_all = 1:7;

% Human-readable subject name for logging.
    subjList = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA', 'CS', 'DT', 'DU', 'RC', 'SR'};
    subjName = subjList{isubj};

% Print run info.
fprintf(' - Subject = %s\n', subjName);
if iscell(namesLocComb) && numel(namesLocComb) >= iLocComb
    locLabel = namesLocComb{iLocComb};
else
    locLabel = 'NA';
end
fprintf(' - Location combo = L%d [%s]\n', iLocComb, locLabel);
fprintf(' - nIter = %d | nJob = %d | iJob = %d\n', nIter, nJob, iJob);
fprintf(' - Fit ModelA = %s\n', strjoin(string(iModelA_fit_all), ' '));
fprintf(' - Fit ModelB = %s\n', strjoin(string(iModelB_fit_all), ' '));
%% Part 1: compute DVs/templates (A1 + A2)
OOD_NOM_Trialwise_compDV_A12(isubj, iLocComb, nIter, nJob, iJob)

%% Part 2: fit NOM for all A/B combinations
for iModelA_fit = iModelA_fit_all
    for iModelB_fit = iModelB_fit_all
        OOD_NOM_Trialwise_fitNOM(isubj, iLocComb, iModelA_fit, iModelB_fit, nIter, nJob, iJob)
    end
end

% End timing.
time_end = datetime('now');
fprintf('%s: OOD_NOM_Human done.\n\n', time_end)
elapsed = time_end - time_start;
fprintf('Time used: %s\n\n\n\n', char(elapsed));

end
