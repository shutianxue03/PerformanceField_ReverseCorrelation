function OOD_NOM_Trialwise_fitNOM(nBasisORI, nBasisSF, flag_whitenDV, isubj, iLocComb, flag_incluCrit, C_contribution, iModelA_fit, iModelB_fit, nIter, nJob, iJob)
%==========================================================================%
% OOD_NOM_Trialwise_fitNOM.m
%--------------------------------------------------------------------------
% Author: Shutian Xue
% Last update: 2025-08-27
%
% Part of the OOD_NOM_Trialwise pipeline.
%
% This function:
% 1) Fits a trial-wise noisy observer model (NOM; Model B) to internal
% variables (IVs) for the TEST set, precomputed in
% OOD_NOM_Trialwise_compIV.m.
% 2) Uses the fitted parameters to predict behavioral metrics (pYES, pC,
% pA) based on binned empirical IVs.
%
% Saved variables (appended to *_compIV.mat):
%       params_est_allIter : [nIter x nParams] fitted parameters
%       nLL_allIter : [nIter x 1] negative log-likelihood
%       pred_metrics_allIter : {nIter x 1} predictions from PR_pred_v7
%
%==========================================================================%

clc; close all;
warning off; % (You may want to remove this once things are stable.)
format compact;

time_start = datetime('now');
fprintf('\n=======================================\n')
fprintf('Part 2: Fit NOM and predict metrics')
fprintf('\n=======================================\n')

% fprintf('%s: Step 2 started.\n\n', time_start)

addpath(genpath('fxn_exp'));
addpath(genpath('fxn_NOM'));
addpath(genpath('fxn_RCplot'));
addpath(genpath('fxn_analysis_RC_v2'));
addpath(genpath('SX_toolbox/bads-master'));

%% General settings (from master config)
%--------------%
SX_RC1_setting; % defines nORI, nSF, namesLocComb, namesModelA, namesModelB, nBins, etc.
%--------------%
flag_fittingStep = 1; % one vs. two step fitting (one step is more standard)
flag_fminconORbads = 2; % 1 = use fmincon (faster, local); 2 = use BADS (slower, more robust)
flag_plot_allIter = 1; % 1 = make summary plots across iterations
flag_plot_perIter = 0; % 1 = plot per-iteration fits (can be slow)
if strcmp('HPC', str_envir), flag_plot_allIter = 0; end % don't plot when running on HPC

%% -------------------- Deterministic RNG (grand seed + per-iteration substreams) -------------------- %%
S_seed = GetGrandSeed(nIter, iJob, nJob, nameFolder_Data);

% One RNG stream for the whole job; each iteration uses its own Substream
stream = RandStream('Threefry', 'Seed', S_seed.grandSeed);
RandStream.setGlobalStream(stream);

% fprintf('%s: seed determined.\n', datetime('now'))

%% Set up file paths and names
if isnumeric(isubj)
    % Human subjects
    subjList = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA', 'CS', 'DT', 'DU', 'RC', 'SR'};
    nblocks_allSubj = [200, 240, 210, 240, 220, 220, 215, 220, 205, 210, 205, 205, 205, 205, 195];

    subjName = subjList{isubj};
    nblocks = nblocks_allSubj(isubj);

    nameFolder_OOD_load = sprintf('%s/%s', nameFolder_Data_OOD, subjName); % behav & energy root
    nameFolder_NOM_save = sprintf('%s/%s/L%d', nameFolder_Data_NOM_Trialwise, subjName, iLocComb);
    nameFolder_Figures_perSubj = sprintf('%s/Human/%s', nameFolder_Figures, subjName);
else
    % Ideal observer or simulated observer
    subjName = isubj{1};
    criterion_DV_true = isubj{2};
    nblocks = 0;

    nameFolder_OOD_load = sprintf('%s/%s', nameFolder_Data_OOD, subjName); % behav & energy root
    nameFolder_NOM_save = sprintf('%s/%s', nameFolder_Data_NOM_Trialwise, subjName);
    nameFolder_Figures_perSubj = sprintf('%s/IO/%s', nameFolder_Figures, subjName);
end

if isempty(dir(nameFolder_NOM_save))
    mkdir(nameFolder_NOM_save);
end

if ~strcmp(str_envir, 'HPC') && isempty(dir(nameFolder_Figures_perSubj))
    mkdir(nameFolder_Figures_perSubj);
end

% File names:
nameFile_compIV = sprintf('%s/n%d_J%d_A%d_compIV', nameFolder_NOM_save, nIter, iJob, iModelA_fit);
nameFile_fitNOM = sprintf('%s/n%d_J%d_A%dB%d', nameFolder_NOM_save, nIter, iJob, iModelA_fit, iModelB_fit);

%% Print header
fprintf(['\nSubject/IO name: %s ' ...
    '\n - L%d [%s]', ...
    '\n - A%d [%s]', ...
    '\n - B%d [%s]', ...
    '\n - Number of iterations = %d\n\n'], ...
    subjName, ...
    iLocComb, namesLocComb{iLocComb}, ...
    iModelA_fit, namesModelA{iModelA_fit}, ...
    iModelB_fit, namesModelB{iModelB_fit}, ...
    nIter);

%% Load trial-wise data and criterion from xx_compIV.mat file
load(nameFile_compIV, 'data_*allIter');
fprintf('%s: Loaded xx_compIV.mat for "data_allIter" .\n\n', datetime('now'))

%% Parameter vectors per Model B
switch flag_incluCrit
    case 1
        % Use the min and max of DV to constrain criterion
        NOMc_lb = min(data_train_allIter{1}.IV);
        NOMc_ub = max(data_train_allIter{1}.IV);
        NOMc_0 = mean([NOMc_lb, NOMc_ub]);
        switch iModelB_fit
            case 1  % FullModel: multi. + additive + shared, criterion
                params0   = [NOMp1_0,   NOMp2_0,   NOMp3_0, NOMc_0];
                params_lb = [NOMp1_lb, NOMp2_lb, NOMp3_lb, NOMc_lb];
                params_ub = [NOMp1_ub, NOMp2_ub, NOMp3_ub, NOMc_ub];

            case 2  % No Nmul: additive + shared (no multiplicative noise)
                params0   = [NOMp2_0,   NOMp3_0, NOMc_0];
                params_lb = [NOMp2_lb, NOMp3_lb, NOMc_lb];
                params_ub = [NOMp2_ub, NOMp3_ub, NOMc_ub];

            case 3  % No Nadd: multi. + shared (no additive noise)
                params0   = [NOMp1_0,   NOMp3_0, NOMc_0];
                params_lb = [NOMp1_lb, NOMp3_lb, NOMc_lb];
                params_ub = [NOMp1_ub, NOMp3_ub, NOMc_ub];

            case 4  % No Nshared: multi. + additive (no shared noise)
                params0   = [NOMp1_0,   NOMp2_0, NOMc_0];
                params_lb = [NOMp1_lb, NOMp2_lb, NOMc_lb];
                params_ub = [NOMp1_ub, NOMp2_ub, NOMc_ub];

            case 5  % Just Nmul: multi. only
                params0   = [NOMp1_0, NOMc_0];
                params_lb = [NOMp1_lb, NOMc_lb];
                params_ub = [NOMp1_ub, NOMc_ub];

            case 6  % Just Nadd: additive only
                params0   = [NOMp2_0, NOMc_0];
                params_lb = [NOMp2_lb, NOMc_lb];
                params_ub = [NOMp2_ub, NOMc_ub];

            case 7  % Just Nshared: shared only
                params0   = [NOMp3_0, NOMc_0];
                params_lb = [NOMp3_lb, NOMc_lb];
                params_ub = [NOMp3_ub, NOMc_ub];

            case 8  % Just criterion: ideal observer (no internal noise)
                params0   = [NOMc_0];
                params_lb = [NOMc_lb];
                params_ub = [NOMc_ub];

            otherwise
                error('OOD_NOM_Trialwise_fitNOM: Unknown iModelB_fit = %d', iModelB_fit);
        end

end

%% Preallocate outputs
nParams = length(namesModelBparams{iModelB_fit});
params_est_allIter = nan(nIter, nParams);
nLL_train_allIter = nan(nIter, 1);
nLL_test_allIter = nLL_train_allIter;
pred_metrics_allIter = cell(nIter, 1);

nameFile_progress = [nameFile_fitNOM, '.mat'];

%% Main estimation loop across iterations
fprintf('%s: Started running %d iterations.\n\n', datetime('now'), nIter)

for iIter = 1:nIter

    fprintf('%d... ', iIter);

    % Deterministic randomness for THIS iteration (global index across jobs)
    stream.Substream = S_seed.iterIdxList(iIter);

    % Extract the training and test set
    data_train = data_train_allIter{iIter}; % (for estimating params)
    data_test = data_test_allIter{iIter}; % for predicting metrics and calculating nLL

    %% Objective function for optimizer: nLL from trial-wise pYES + pairwise pA
    switch flag_fittingStep

        case 1 % one-step fitting (more standard)
            iStep = 0; % iStep=0: use the full joint nLL (pYES + pA)
            switch flag_incluCrit
                % case 0
                %     %------------------------------%
                %     fxn_estParams = @(paramsNOM) fxn_getError_v7(iModelB_fit, paramsNOM, data_train, criterion_z, iStep); % the full model has multi., additive noise (shared and additive across passes)
                %     %------------------------------%
                case 1
                    %------------------------------%
                    fxn_estParams = @(paramsNOM) fxn_getError_v8(iModelB_fit, paramsNOM, data_train, C_contribution, iStep); % the full model has multi., additive noise (shared and additive across passes)
                    %------------------------------%
            end

            % Estimate parameters
            if flag_fminconORbads == 1 % Faster, local search
                [params_est, nLL_train] = fmincon(fxn_estParams, params0, [], [], [], [], params_lb, params_ub, [], options_fmin);
            else % BADS: more robust global + local search
                [params_est, nLL_train] = bads(fxn_estParams, params0, params_lb, params_ub, [], [], [], options_bads);
            end

        case 2 % Two-stage fitting: stabilize models that include a shared component across two passes
            if any(iModelB_fit == [1,2,3,7])

                % ==== Step 1: fit the "base" noise terms (excluding the shared term) ====
                iStep = 1;
                switch flag_incluCrit
                    % case 0
                    %     %------------------------------%
                    %     fxn_estParams = @(paramsNOM) fxn_getError_v7(iModelB_fit, paramsNOM, data_train, criterion_z, iStep);
                    %     %------------------------------%
                    case 1
                        %------------------------------%
                        fxn_estParams = @(paramsNOM) fxn_getError_v8(iModelB_fit, paramsNOM, data_train, C_contribution, iStep);
                        % NOTE: this line is just building the fitting fxn, so no need to preprocess paramsNOM here
                        %------------------------------%
                end

                % Index the parameters to be included in step 1 (calculating nLL for detection prob.)
                switch flag_incluCrit
                    % case 0
                    %     switch iModelB_fit
                    %         case 1, indStep1=1:2;
                    %         case 2, indStep1=1;
                    %         case 3, indStep1=1;
                    %         case 4, indStep1=1;
                    %     end
                    case 1
                        switch iModelB_fit
                            case 1, indStep1=[1,2,4]; % params0=[Nmul,Nadd,Nshared,criterion] -> pick [Nmul,Nadd,criterion]
                            case 2, indStep1=[1,3]; % params0=[Nadd,Nshared,criterion] -> pick [Nadd,criterion]
                            case 3, indStep1=[1,3]; % params0=[Nmul,Nshared,criterion] -> pick [Nmul,criterion]
                            case 4, indStep1=[1,3]; % params0=[Nmul,Nadd,criterion] -> pick [Nmul,criterion] (unused in two-step)
                        end
                end

                % Adjust parameter range
                params0_step1 = params0(indStep1);
                params_lb_step1 = params_lb(indStep1);
                params_ub_step1 = params_ub(indStep1);

                % Estimate parameters
                if flag_fminconORbads == 1 % Faster, local search
                    [params_est_fromStep1, nLL_step1] = fmincon(fxn_estParams, params0_step1, [], [], [], [], params_lb_step1, params_ub_step1, [], options_fmin);
                else % BADS: more robust global + local search
                    [params_est_fromStep1, nLL_step1] = bads(fxn_estParams, params0_step1, params_lb_step1, params_ub_step1, [], [], [], options_bads);
                end

                % ==== Step 2. Fix the Step-1 estimates and fit the shared term ====
                iStep = 2;
                switch flag_incluCrit
                    % case 0
                    %     %------------------------------%
                    %     fxn_estParams = @(paramsNOM) fxn_getError_v7(iModelB_fit, paramsNOM, data_train, criterion_z, iStep, params_est_fromStep1); % the full model has multi., additive noise (shared and additive across passes)
                    %     %------------------------------%
                    case 1
                        %------------------------------%
                        fxn_estParams = @(paramsNOM) fxn_getError_v8(iModelB_fit, paramsNOM, data_train, C_contribution, iStep, params_est_fromStep1); % the full model has multi., additive noise (shared and additive across passes)
                        %------------------------------%
                end

                % Index the parameters to be included in step 2 (calculating nLL for detection prob.)
                switch flag_incluCrit
                    % case 0
                    %     switch iModelB_fit
                    %         case 1, indStep2=3;
                    %         case 3, indStep2=2; % B3=[sigmaAdd,sigmaShared]; step2 fits sigmaShared (idx 2)
                    %         case 4, indStep2=2; % B4=[Nmul,sigmaShared]; step2 fits sigmaShared (idx 2)
                    %     end
                    case 1
                        switch iModelB_fit
                            case 1, indStep2=3; % params0=[Nmul,Nadd,Nshared,criterion] -> pick Nshared
                            case 2, indStep2=2; % params0=[Nadd,Nshared,criterion] -> pick Nshared
                            case 3, indStep2=2; % params0=[Nmul,Nshared,criterion] -> pick Nshared
                                % case 4, indStep2=nan; % iModelB=4 has no shared term, so not in two-step set
                        end
                end

                % Adjust parameter range
                params0_step2 = params0(indStep2);
                params_lb_step2 = params_lb(indStep2);
                params_ub_step2 = params_ub(indStep2);

                % Estimate parameters
                if flag_fminconORbads == 1 % Faster, local search
                    [params_est_step2, nLL_step2] = fmincon(fxn_estParams, params0_step2, [], [], [], [], params_lb_step2, params_ub_step2, [], options_fmin);
                else % BADS: more robust global + local search
                    [params_est_step2, nLL_step2] = bads(fxn_estParams, params0_step2, params_lb_step2, params_ub_step2, [], [], [], options_bads);
                end

                % Combine est. params of two steps
                switch flag_incluCrit
                    % case 0
                    %     switch iModelB_fit
                    %         case 1, params_est = [params_est_fromStep1([1,2]), params_est_step2]; % B1: [Nmul,Nadd,Nshared]
                    %         case 3, params_est = [params_est_fromStep1(1), params_est_step2]; % B3: [sigmaAdd,sigmaShared]
                    %         case 4, params_est = [params_est_fromStep1(1), params_est_step2]; % B4: [Nmul,sigmaShared]
                    %     end
                    case 1
                        switch iModelB_fit
                            case 1, params_est = [params_est_fromStep1([1,2]), params_est_step2, params_est_fromStep1(3)]; % [Nmul, Nadd, Nshared, criterion]
                            case 2, params_est = [params_est_fromStep1(1), params_est_step2, params_est_fromStep1(2)]; % [Nadd, Nshared, criterion]
                            case 3, params_est = [params_est_fromStep1(1), params_est_step2, params_est_fromStep1(2)]; % [Nmul, Nshared, criterion]
                        end
                end
                nLL_train = nLL_step1 + nLL_step2;

            else % For models without shared term (e.g., B=4), fit all params together
                iStep = 0;
                switch flag_incluCrit
                    % case 0
                    %     %------------------------------%
                    %     fxn_estParams = @(paramsNOM) fxn_getError_v7(iModelB_fit, paramsNOM, data_train, criterion_z, iStep); % the full model has multi., additive noise (shared and additive across passes)
                    %     %------------------------------%
                    case 1
                        %------------------------------%
                        fxn_estParams = @(paramsNOM) fxn_getError_v8(iModelB_fit, paramsNOM, data_train, C_contribution, iStep); % the full model has multi., additive noise (shared and additive across passes)
                        %------------------------------%
                end
                % Estimate parameters
                if flag_fminconORbads == 1 % Faster, local search
                    [params_est, nLL_train] = fmincon(fxn_estParams, params0, [], [], [], [], params_lb, params_ub, [], options_fmin);
                else % BADS: more robust global + local search
                    [params_est, nLL_train] = bads(fxn_estParams, params0, params_lb, params_ub, [], [], [], options_bads);
                end
            end
    end


    %% Caluclate nLL and predict metrics
    iStep = 0; % !iStep=0!: use the full joint nLL (pYES + pA).
    switch flag_incluCrit
        % case 0
        %     %------------------------------%
        %     % pred_train = PR_pred_v7(iModelB_fit, params_est, data_train, criterion_z, nBins, flag_plot_perIter);
        %     %------------------------------%

        %     % Calculate nLL for the test set
        %     %------------------------------%
        %     nLL_test = fxn_getError_v7(iModelB_fit, params_est, data_test, criterion_z, iStep);
        %     %------------------------------%

        %     % Predict binned metrics from estimated parameters
        %     %------------------------------%
        %     pred_test = PR_pred_v7(iModelB_fit, params_est, data_test, criterion_z, nBins, flag_plot_perIter);
        %     %------------------------------%

        case 1
            %------------------------------%
            % pred_train = PR_pred_v8(iModelB_fit, params_est, data_train, nBins, flag_plot_perIter);
            %------------------------------%

            % Calculate nLL for the test set
            %------------------------------%
            nLL_test = fxn_getError_v8(iModelB_fit, params_est, data_test, C_contribution, iStep);
            %------------------------------%

            % Predict binned metrics from estimated parameters
            %------------------------------%
            pred_test = PR_pred_v8(iModelB_fit, params_est, data_test, nBins, flag_plot_perIter);
            %------------------------------%
    end

    %% Compile
    nLL_train_allIter(iIter) = nLL_train;
    nLL_test_allIter(iIter) = nLL_test;
    params_est_allIter(iIter, :) = params_est(:).';
    pred_metrics_allIter{iIter} = pred_test;

    %% Save a progress report in the folder to indicate the finished iteration and time spent
    time_progress = ceil(minutes(datetime('now')-time_start)); % round up to minutes
    save(nameFile_progress, 'iIter')
    % rename (to avoid saving one file for each iteration)
    nameFile_progress_new = sprintf('%s_%d_%dmin.mat', nameFile_fitNOM, iIter, time_progress);
    movefile(nameFile_progress, nameFile_progress_new);
    nameFile_progress = nameFile_progress_new;

end % end of iIter

fprintf('\n\n%s: All iterations done.\n\n', datetime('now'))

%% Save results (append onto *_compIV.mat)
save(nameFile_fitNOM, '*_allIter', 'time_progress');

% Delete the progress report
if exist(nameFile_progress_new, 'file')
    delete(nameFile_progress_new);
end

fprintf('%s: Outputs saved.\n\n', datetime('now'))

%% Plot
if flag_plot_allIter
    %-------------------%
    NOMplot_fitNOM;
    %-------------------%
    fprintf('\n\n%s: Plots created.\n\n', datetime('now'))
end
close all;

%% Timing info
time_end = datetime('now');
fprintf('%s: Fit NOM done.\n\n', time_end)
elapsed = time_end - time_start;
fprintf('Time used: %s\n\n\n\n', char(elapsed));


end
