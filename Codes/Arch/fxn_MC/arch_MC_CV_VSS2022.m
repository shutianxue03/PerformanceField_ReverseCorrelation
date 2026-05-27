% 
% % 
% clc
% 
% ifeature = 2;
% ifamily = 2;
% isubj = 1;
% ncv = 2; % number of rounds of cross validation
% nrep = 2; % number of reps to find the global minimum
% 
% %%
% addpath(genpath('Data_MC'))
% addpath(genpath('Data_OOD'))
% addpath(genpath('fxn_analysis_RC'))
% 
% load('Data_OOD/n10.mat')
% colors2Loc = [ 1 0 0; 0 0 1];
% 
% k = 10; % k-fold cross validation
% subjName = subjList{isubj};
% 
% fprintf('Model comparison starts:\n   [%s - M%d] %s\n', namesFeature{ifeature}, ifamily, namesFamily_all{ifamily})
% fprintf('   Subj#%d %s\n   %d rounds in the %d-fold cross validation\n ', isubj, subjName, ncv, k)
% fprintf('   %d reps to find the global minimum\n\n', nrep)
% 
% %%
% sampParam0 = @(lb, ub, n) rand(1,n) .* (ub-lb) + lb;
% options = optimoptions('fmincon','MaxIterations',5000,'Display','off');
% nitp = 2; % number of multiples of interpolated pints
% nLocTrue = 2;
% iLocTrue = [1,8];
% 
% %% extract data
% % extract y data
% yData = squeeze(margK_allSubj{ifeature}(isubj, 3, 1, :, :));
% 
% % interpolate x-axis
% if ifeature==2, xaxis_inUse = 2.^xaxis_all{ifeature};
% else, xaxis_inUse = xaxis_all{ifeature};
% end
% nfilters = length(xaxis_inUse);
% xaxis_itp = interp1(1:nfilters, xaxis_inUse, linspace(1, nfilters, nfilters*nitp),'spline');
% nx = nfilters*nitp;
% 
% %% model setting
% % names of params
% namesParams = namesParams_all{ifamily};
% nparams_full = length(namesParams);
% 
% % decide number of models & the param index (1 = specified, 0 = shared)
% if (ifamily == 3) || (ifamily == 7) % the truncation term is always left free
%     nmodels = 2^(nparams_full-1);
%     paramInd_all_ = fxn_getParamInd(nparams_full-1);
%     paramInd_all = [paramInd_all_, ones(nmodels,1)];
% else
%     nmodels = 2^nparams_full;
%     paramInd_all = fxn_getParamInd(nparams_full);
% end
% 
% % extract ub and lb
% ub_full = ub_full_all{ifamily};
% lb_full = lb_full_all{ifamily};
% 
% %% ETC
% disp(datetime('now'))
% ETC = round(3.75*ncv*nrep*nmodels/60);
% fprintf('ETC: %.1f minutes\n', ETC);
% 
% %% cross-validation setting
% hpartition = cvpartition(nx, 'kfold', k);
% TestSize = hpartition.TestSize;
% testSize_cumsum = cumsum(TestSize);
% 
% %%
% RSS_ML_ave_all = nan(ncv, nmodels);
% RSS_ML_std_all = nan(ncv, nmodels);
% 
% %%
% parfor icv = 1:ncv
%     rng shuffle
%     ind = [randperm(nx, nx); randperm(nx, nx)]; % decide the random sequence of this round
%     
%     RSS_LS_all = nan(k, nmodels);
%     RSS_ML_all = RSS_LS_all;
%     
%     for ik = 1:k
%         if ik==1, istart = 1; else, istart = TestSize(ik-1)+1;end
%         iend = testSize_cumsum(ik);
%         
%         ntest = length(istart:iend);
%         ntrain = nx - ntest;
%         
%         % empty containers
%         x_test = nan(2, ntest);
%         y_test  = x_test;
% 
%         x_train = nan(2,ntrain);
%         y_train = x_train;
% 
%         for iline = 1:2
%             % get the index for test group
%             ind_test = ind(iline, istart:iend);
%             ind_train = 1:nx;
%             ind_train(ind_test) = [];
%             % x
%             x_test(iline, :) = xaxis_itp(ind_test);
%             x_train(iline, :) = xaxis_itp(ind_train);
%             % y
%             y_test(iline, :) = yData(iline, ind_test);
%             y_train(iline, :) = yData(iline, ind_train);
%             
%         end % end of iline
%         
% %         figure('Position', [3000 600 1800 1300])
%         for imodel = 1:nmodels % loop thru param combinations
%             
%             paramInd = paramInd_all(imodel, :);
%             nparamsAll = sum(nLocTrue.^paramInd);
%             [ub, lb, params0] = MC_getLimits(paramInd, ub_full, lb_full, nLocTrue);
%             % weight error terms by trial number
%             %                 if nsubj>1, ntrialsProp = ntrials_allSubj'/sum(ntrials_allSubj); else, ntrialsProp = 1; end
%             
%             % fit model using LS method
% %             fitMode = 1;
% %             fxn_getDev_LS = @(kernelParams) MC_getDev_CV(paramInd, kernelParams, x_train, y_train, ifamily, fitMode);
% %             problem_LS = createOptimProblem('fmincon','objective', fxn_getDev_LS,'x0',params0,'lb',lb,'ub',ub,'options',options);
% %             ms_LS = MultiStart('StartPointsToRun', 'bounds', 'Display', 'off');
% %             [params_est_LS] = fxn_parfor(ms_LS, problem_LS, nrep);
% %             % validate
% %             RSS_LS_all(ik, imodel) = MC_getDev_CV(paramInd, params_est_LS, x_test, y_test, ifamily, fitMode);
%             
%             % ML method
%             fitMode = 2;
%             % calibrate
%             fxn_getDev_ML = @(kernelParams) MC_getDev_CV(paramInd, kernelParams, x_train, y_train, ifamily, fitMode);
%             problem_ML = createOptimProblem('fmincon','objective', fxn_getDev_ML,'x0',params0,'lb',lb,'ub',ub,'options',options);
%             ms_ML = MultiStart('StartPointsToRun', 'bounds', 'Display', 'off');
%             [params_est_ML] = fxn_parfor(ms_ML, problem_ML, nrep);
%             %  validate
%             RSS_ML_all(ik, imodel) = MC_getDev_CV(paramInd, params_est_ML, x_test, y_test, ifamily, fitMode);
%             
%             % plot the predicted value
%             params_est = params_est_ML;
%             y_testpred = RC_predKernel_MC(paramInd, x_test, 2, nparams_full, params_est, ifamily);
%             yData_pred = RC_predKernel_MC(paramInd, xaxis_itp, 2, nparams_full, params_est, ifamily);
% %             MCplot_temp
%         end % end of imodel (parfor)
%        
%     end % end of ik
%     RSS_ML_ave_all(icv, :) = mean(RSS_ML_all,1);
%     RSS_ML_std_all(icv, :) = std(RSS_ML_all, [], 1);
%     fprintf('Finished icv = %d\n', icv)
% end % end of icv (cross validation)
% 
% disp(datetime('now'))
% 
% %% some calculation
% [~, ibest] = min(mean(RSS_ML_ave_all));
% fprintf('%s: best model is #%d, param Ind: %s\n', subjName, ibest, num2str(paramInd_all(ibest, :)))
% 
% %% save
% save(sprintf('Data_MC/CV_%s_%s_M%d', subjList{isubj}, namesFeature{ifeature}, ifamily), 'RSS*', 'ibest')
% 
% %%
% function [params_est_LS, RSS, flag, output, solution] = fxn_parfor(ms_LS, problem_LS, nrep)
%  [params_est_LS, RSS, flag, output, solution]= run(ms_LS, problem_LS, nrep);
% end
% 
