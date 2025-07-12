% % estimate internal-external noise ratio by minimizing the calculated pA and estimated pA
% % only analyze data, not plotting
% 
% clc
% warning off
% format compact
% 
% addpath(genpath('Data_OOD'))
% addpath(genpath('fxn_model'))
% addpath(genpath('fxn_analysis_RC'))
% 
% %% shell
% isubj = 3; % the selected subj among {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA', 'CS', 'DT'};
% iLoc_all = [1,8];
% templateType.i = 3;
% convolveType.i = 1;
% IVType.i = 1;
% normIV = 1;
% errorComp = 'dprime, criterion, pA3, pC3';
% 
% ni = 100; % number of bootstrapping, 100
% nrep = 10; % 20
% nparams_model = 3; % 1 = alpha (assumed the same for PRS and ABS)
%                                % 2 = alpha (assumed the same for PRS and ABS) and thresh
%                                % 3 = alpha_PRS, alpha_ABS and thresh
% errorFxn = '@(pred, data) sum(-log(normpdf(pred, data, 1)), ''all'')';
% 
% 
% %% PARAMS
% subjList = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA', 'CS', 'DT'};
% nblocks_allSubj = [200, 240, 210, 240, 220, 220, 215, 220, 205, 210, 205, 205];
% 
% subjName = subjList{isubj};
% nblocks = nblocks_allSubj(isubj);
% 
% nB = 1;
% ntrials_sim = 1e5;
% nfiltersOri = 29;
% nfiltersSF = 29;
% itype = 3;
% nLoc = length(iLoc_all);
% options = optimoptions('fmincon','MaxIterations', 1e4, 'Display','off');
% fprintf('\nNumber of boostraps (nB): %d\nNumber of ORI filters: %d\nNumber of SF filters: %d\n\n', nB, nfiltersOri, nfiltersSF)
% 
% fileName_perSubj = sprintf('Data_OOD/%s_model.mat', subjName);
% fileDir_perSubj = dir(fileName_perSubj);
% 
% fprintf('\nSubj name: %s \nNumber of blocks (nblocks): %d\n', subjName, nblocks)
% 
% %%  define the range
% % Params:
% % 1. internal-external noise ratio of the PRS trials
% % 2. internal-external noise ratio of the ABS trials
% % 3. criterion (not the criterion between -1 to 1)
% 
% % Define the initial/lb/ub
% switch nparams_model
%     case 1
%         params0 = .5;
%         params_lb = 0;
%         params_ub = 2;
%     case 2
%         params0 = [.5, 0];
%         params_lb = [0, -1];
%         params_ub = [2, 1];
%     case 3
%         params0 = [.5, .5, 0];
%         params_lb = [0, 0, -1];
%         params_ub = [2, 1.5, 1];
% end
% assert(nparams_model == length(params_ub));
% 
% params0 = repmat(params0, [nLoc,1]);
% params_lb = repmat(params_lb, [nLoc,1]);
% params_ub = repmat(params_ub, [nLoc,1]);
% load('Data_OOD/filtersFile', 'energy2D_gabor')
% 
% %% save the description of the analysis
% templateType.all = {'(1) the 2D kernel (raw)',...
%                              '(2) the 2D kernel (all positive)', ...
%                              '(3) the reconstructed 2D kernel (all positive)', ...
%                              '(4) the energy profile of the gabor (CST=1)'};
% convolveType.all = {'(1) The 2D template is convolved with the 2D energy profile'};
% IVType.all = {'(1) sum all channels up', ...
%                     '(2) select the highest value'};
% 
% README = sprintf('\n Loc fitted = %s\n ni = %d\n nrep = %d\n Template: %s\n Convolve type = %s\n Sum up or select the max = %s\n Normalize IV distributions (1=YES, 0=NO): %d\n Criterion [initial/lb/ub]: %.1f, %.1f, %.1f\n Error fxn: %s\n Error components: %s\n \n', ...
%     num2str(iLoc_all), ni, nrep, ...
%     templateType.all{templateType.i}, ...
%     convolveType.all{convolveType.i}, ...
%     IVType.all{IVType.i}, ...
%     normIV, ... % whether IV distributions are normalized?
%     [params0(1,3), params_lb(1,3), params_ub(1,3)], ...
%     errorFxn,...
%     errorComp); % error components
% fprintf(README)
% 
% subjName_ = subjName;
% 
% %%
% if isempty(fileDir_perSubj)
%     
%     load(sprintf('Data_OOD/%s%d_B%d_%d_%d.mat', subjName_, nblocks, nB, nfiltersOri, nfiltersSF))
%     load(sprintf('Data_OOD/%s%d_energy_%d_%d.mat', subjName_, nblocks, nfiltersOri, nfiltersSF))
%     load(sprintf('Data_OOD/%s%d_behavMeas.mat', subjName_, nblocks))
%     
%     nAllTrials = length(energy2D_allT_perComb{1});
%     ntrials = nAllTrials/2;
%     
%     nsubj = 1; % do NOT change
%     subjName = subjName_;
%     
%     %% 1. calculate the internal variable
%     fprintf('1. Calculating the internal variable... ')
%     
%     %%%%%%%
%     fxn_getIV
%     %%%%%%%
%     
%     % compile data
%     data.IV_PRS = IV_PRS_perLoc;
%     data.IV_ABS = IV_ABS_perLoc;
%     data.dprime = dprime_allT_perComb(iLoc_all);
%     data.criterion = criterion_allT_perComb(iLoc_all);
%     data.pC = pC_allT_perComb(iLoc_all, 3);
%     data.pHit = pC_allT_perComb(iLoc_all, 1);
%     data.pFA = pC_allT_perComb(iLoc_all, 2);
%     data.pA = pA_allT_perComb(iLoc_all, 3);
%     data.pA_PRS = pA_allT_perComb(iLoc_all, 1);
%     data.pA_ABS = pA_allT_perComb(iLoc_all, 2);
%     
%     fprintf('Done!\n')
%     
%     %% 2. fit & predict
%     fprintf('2. Fitting...')
%     
%     params_est_all = nan(ni, nLoc, nparams_model);
%     dprime_pred_all = nan(ni, nLoc);
%     criterion_pred_all = dprime_pred_all;
%     pA_pred_all = dprime_pred_all;
%     pA_PRS_pred_all = dprime_pred_all;
%     pA_ABS_pred_all = dprime_pred_all;
%     pC_pred_all = dprime_pred_all;
%     pHit_pred_all = dprime_pred_all;
%     pFA_pred_all = dprime_pred_all;
%     
%     tStart = tic;
%     for ii = 1:ni
%         
%         %%%%%%%
%         %   Fitting  %
%         %%%%%%%
%         fxn_estParams = @(params) fxn_getError(params, data, ntrials_sim, normIV);
%         problem_ML = createOptimProblem('fmincon','objective', fxn_estParams,'x0', params0, 'lb', params_lb, 'ub', params_ub, 'options',options);
%         ms_ML = MultiStart('StartPointsToRun', 'bounds', 'Display', 'off', 'UseParallel', 1);
%         
%         params_est = run(ms_ML, problem_ML, nrep);
% %         params_est = fmincon(fxn_estParams, params0, [], [], [], [], params_lb, params_ub, [],options);
%         params_est_all(ii, :, :) = params_est;
%         
%         %%%%%%%%%
%         %   Prediction  %
%         %%%%%%%%%
%         pred = fxn_getEstimates(params_est, data.IV_PRS, data.IV_ABS, ntrials_sim, normIV, 0);
%         dprime_pred_all(ii, :) = pred.dprime;
%         criterion_pred_all(ii, :) = pred.criterion;
%         pC_pred_all(ii, :) = pred.pC;
%         pHit_pred_all(ii, :) = pred.pHit;
%         pFA_pred_all(ii, :) = pred.pFA;
%         pA_pred_all(ii, :) = pred.pA;
%         pA_PRS_pred_all(ii, :) = pred.pA_PRS;
%         pA_ABS_pred_all(ii, :) = pred.pA_ABS;
%         
%         if ~mod(ii, ni/10), fprintf('='), end
%     end
%     
%     tEnd = toc(tStart);
%     fprintf('Done! %.1f mins\n', tEnd/60)
%     
%     %% compile
%     pred_all.dprime= dprime_pred_all;
%     pred_all.criterion= criterion_pred_all;
%     pred_all.pC = pC_pred_all;
%     pred_all.pHit = pHit_pred_all;
%     pred_all.pFA = pFA_pred_all;
%     pred_all.pA = pA_pred_all;
%     pred_all.pA_PRS = pA_PRS_pred_all;
%     pred_all.pA_ABS = pA_ABS_pred_all;
%     
%     %% 4. save data
%     save(fileName_perSubj, 'template*', 'README', 'data','*params_est_all', '*pred_all')
%     fprintf(' ========== Finished ==========\n\n')
%     
% else
%     fprintf('File exists.\n')
%     
% end
% 
% 
% 
