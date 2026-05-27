function OOD_NOM_trialWise_arch(isubj, iLocComb, iModelA_fit_all, nIter, nJob, iJob)
% OOD_NOM_Trialwise
% Wrapper that follows the OOD_sim two-stage flow:
%   1) OOD_NOM_Trialwise_compIV_A12
%   2) fit stage via shell_NOM_fitNOM.sh (or direct MATLAB fallback)

clc; close all;
format compact;

if nargin < 1 || isempty(isubj), isubj = 1; end
if nargin < 2 || isempty(iLocComb), iLocComb = 1; end
if nargin < 3 || isempty(iModelA_fit_all), iModelA_fit_all = 1; end
if nargin < 4 || isempty(nIter), nIter = 200; end
if nargin < 5 || isempty(nJob), nJob = 1; end
if nargin < 6 || isempty(iJob), iJob = 1; end

% Keep these aligned with the NOM trialwise pipeline defaults used by wrappers.
nBasisORI = 6;
nBasisSF = 6;
flag_whitenDV = 1;
lambda_whiten = 1;
flag_regressType = 2;

time_start = datetime('now');
fprintf('\n=======================================\n')
fprintf('Run OOD_NOM Trialwise Pipeline')
fprintf('\n=======================================\n')
fprintf('Subject: %s\n', mat2str(isubj))
fprintf('Location: %d\n', iLocComb)
fprintf('A models: %s\n', strjoin(string(iModelA_fit_all), ' '))
fprintf('nIter=%d, iJob=%d/%d\n\n', nIter, iJob, nJob)

% Step 1: compute IV/template files (A1 + A2) using the updated A12 pipeline.
OOD_NOM_Trialwise_compIV_A12( ...
    isubj, iLocComb, lambda_whiten, ...
    flag_whitenDV, nIter, nJob, iJob);

% Step 2: run fit stage through the shell wrapper.
thisFile = mfilename('fullpath');
[scriptDir, ~, ~] = fileparts(thisFile);
shellFit = fullfile(scriptDir, 'shell_NOM_fitNOM.sh');

if exist(shellFit, 'file') == 2 && system('command -v sbatch >/dev/null 2>&1') == 0
    for iModelA_fit = iModelA_fit_all
        cmd = sprintf('sbatch "%s" %d %d %d %d %d %d', ...
            shellFit, isubj, iLocComb, iModelA_fit, nIter, nJob, iJob);
        fprintf('Submitting fit stage: %s\n', cmd);
        [status, out] = system(cmd);
        if status ~= 0
            error('Failed to submit fit stage.\nCommand: %s\nOutput:\n%s', cmd, out);
        end
        fprintf('%s\n', strtrim(out));
    end
else
    % Fallback for non-Slurm environments: run fit directly in MATLAB.
    iModelB_fit_all = 1:7;
    for iModelA_fit = iModelA_fit_all
        for iModelB_fit = iModelB_fit_all
            OOD_NOM_Trialwise_fitNOM( ...
                isubj, iLocComb, flag_whitenDV, iModelA_fit, iModelB_fit, ...
                nIter, nJob, iJob);
        end
    end
end

time_end = datetime('now');
fprintf('\n%s: OOD_NOM_Trialwise wrapper done.\n', time_end)
fprintf('Time used: %s\n\n', char(time_end - time_start));

end% function OOD_NOM_Trialwise(isubj, iLocComb, iModelA, nIterations)
% % OOD_NOM_trialWise: This function generates IVs and fits different models to them.
% % Inputs:
% %      isubj: index of subj
% %      iLocComb: 1=Fovea, 8=periF(6 deg ecc), 6=HM, 7=VM, 5=LVM, 3=UVM
% %      iModelA: 1=core model, 2=permuted template, 3=use IO template
% %      ni: number of iterations

% clc
% rng(123);   % define see for reproducibility
% time_start = datetime('now');

% %% Generate IVs
% fprintf('\n\n=====================\nGenerating IVs (Model A%d)\n=====================\n\n', iModelA)

% % The settings below can be placed in SX_RC1_setting
% % IVType = 1;            % 1=sum of the dot product/convolution; 2=max; 3=normalized
% % templateType = 3; % (1) raw (2) reconstructed kernel (3) mirrored template
% % itype_template = 2; % 1=estimate template from PRS trials, ABS trials, or BOTH trials
% % flag_PatchMode = 1; % if flag_PatchMode == 1, patchMode = 'T'; else, patchMode = 'N'; end

% %------------------------------%
% OOD_NOM_Trialwise_compIV(isubj, iLocComb, iModelA, nIterations);
% %------------------------------%

% %% Fit different models to IVs and make predictions
% iModelB_all = 1:7; % Check the full list of models in SX_RC1_setting, search "namesModelB"
% % iModelB_all = 2;
% for iModelB = iModelB_all
%     fprintf('\n\n=====================\nModel B%d\n=====================\n', iModelB)
%     % ------------------------------%
%     OOD_NOM_Trialwise_fitNOM(isubj, iLocComb, iModelA, iModelB, nIterations)
%     % ------------------------------%
% end

% %% Ending
% time_end = datetime('now');
% elapsed = time_end - time_start;
% fprintf('\n\n*** TOTAL time used: %s)\n', char(elapsed));
