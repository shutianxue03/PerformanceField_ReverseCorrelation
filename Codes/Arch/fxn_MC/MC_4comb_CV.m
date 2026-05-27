% 
% % fit the selected model to ALL observers, weighted by number of trials
% % using cross-validation!!
% 
% 
% 
% ifeature=2;
% ifamily=3; 
% 
% sampParam0 = @(lb, ub, n) rand(1,n) .* (ub-lb) + lb;
% options = optimoptions('fmincon','MaxIterations',5000,'Display','off');
% nrep = 2; % number of reps to find the global minimum
% nitp = 2; % number of multiples of interpolated pints
% ncomb=1;
% ntypes = 1; itype = 3;
% nLocTrue = 2;
% k = 10; % k-fold cross validation
% ncv = 1; % number of rounds of cross validation
% 
% 
% %%
% % names of params
% namesParams = namesParams_all{ifamily};
% nparams_full = length(namesParams);
% % x-axis
% if ifeature==2, xaxis_inUse = 2.^xaxis_all{ifeature};
% else, xaxis_inUse = xaxis_all{ifeature};
% end
% nfilters = length(xaxis_inUse);
% xaxis_itp = interp1(1:nfilters, xaxis_inUse, linspace(1, nfilters, nfilters*nitp),'spline');
% nx = nfilters*nitp;
% ntrain = round(nx/k*(k-1));
% ntest = nx - ntrain;
% hpartition = cvpartition(nx, 'kfold', k);
% testSize_cumsum = cumsum(hpartition.TestSize);
% 
% %%
% % number of models
% if (ifamily == 3) || (ifamily == 7) % the truncation term is always left free
%     nmodels = 2^(nparams_full-1);
% else
%     nmodels = 2^nparams_full;
% end
% 
% % extract ub and lb
% ub_full = ub_full_all{ifamily};
% lb_full = lb_full_all{ifamily};
% 
% % decide the param index (1 = specified, 0 = shared)
% if (ifamily == 3) || (ifamily == 7) % for model with a truncation term, it's always left free to vary
%     paramInd_all_ = fxn_getParamInd(nparams_full-1);
%     paramInd_all = [paramInd_all_, ones(nmodels,1)];
% else
%     paramInd_all = fxn_getParamInd(nparams_full);
% end
% 
% % empty containers
% yData_allSubj_all = cell(ncomb, ntypes);
% yPred_LS_allSubj_all = yData_allSubj_all;
% yPred_ML_allSubj_all = yData_allSubj_all;
% 
% % estimate time of completion
% % ETC = round(12 * ncomb * ntypes * nsubj *nmodels/60/4);
% % fprintf('Model comparison starts\n%s [M%d-%s], nmodels=%d. ETC: %d minutes.\n', namesFeature{ifeature}, ifamily, namesFamily_all{ifamily}, nmodels, ETC)
% disp(datetime('now'))
% 
% for icomb = 1:ncomb
%     
%     % decide the Loc index of the pair to be compared
%     iLocTrue = combInd(icomb, :);
%     
%     for itype = 3%1:ntypes
%         
%         % extract y data
%         yData_allSubj = squeeze(margK_allSubj{ifeature}(:, itype, icomb, :, :));
%         
%         % empty containers
%         yPred_LS_allSubj = nan(nsubj, nmodels, 2, nfilters*nitp);
%         yPred_ML_allSubj = yPred_LS_allSubj;
%         RSS_LS_ave_allSubj = nan(nsubj, ncv, nmodels);
%         RSS_LS_std_allSubj = nan(nsubj, ncv, nmodels);
%         RSS_best_allSubj = nan(1, nsubj);
%         ibest_LS_allSubj = nan(1, nsubj);
%         
%         % fit each subj
%         for isubj = 1:nsubj
%             
%             % extract data
%             yData = squeeze(yData_allSubj(isubj, :, :));
%             
%             for icv = 1:ncv
%                 ind = [randperm(nx, nx);randperm(nx, nx)];
%                 RSS_LS_ave_all = nan(ncv, nmodels);
%                 RSS_LS_std_all = nan(ncv, nmodels);
%                 for ik = 1:k
%                     if ik==1, istart = 1; else, istart = hpartition.TestSize(ik-1)+1;end
%                     iend = testSize_cumsum(ik);
%                     
%                     % empty containers
%                     ntest = length(istart:iend);
%                     ntrain = nx - ntest;
%                     yData_train = nan(2,ntrain);
%                     x_train = yData_train;
%                     ind_test= nan(2, ntest);
%                     x_test = ind_test;
%                     yData_test  = ind_test;
%                     
%                     for iline = 1:2
%                         % get the index for test group
%                         ind_test(iline, :) = ind(iline, istart:iend);
%                         % x
%                         x_test(iline, :) = xaxis_itp(ind_test(iline, :));
%                         x_train_ = xaxis_itp;
%                         x_train_(ind_test(iline, :)) = [];
%                         x_train(iline, :) = x_train_;
%                         % y
%                         yData_test(iline, :) = yData(iline, ind_test(iline, :));
%                         yData_train_ = yData(iline, :);
%                         yData_train_(ind_test(iline, :)) = [];
%                         yData_train(iline, :) = yData_train_; %Get rid of the element supposed to be deleted
%                     end % end of iline
%                     
%                     %             %  empty containers
%                     %             yPred_LS_allModels = nan(nmodels, 2, nfilters*nitp);
%                     %             yPred_ML_allModels = nan(nmodels, 2, nfilters*nitp);
%                     
%                     for imodel = 1:nmodels % loop thru param combinations
%                         
%                         paramInd = paramInd_all(imodel, :);
%                         nparamsAll = sum(nLocTrue.^paramInd); 
%                         [ub, lb, params0] = MC_getLimits(paramInd_all(imodel, :), ub_full, lb_full, nLocTrue);
%                         % weight error terms by trial number
%                         %                 if nsubj>1, ntrialsProp = ntrials_allSubj'/sum(ntrials_allSubj); else, ntrialsProp = 1; end
%                         
%                         % fit model ML method (minimizing RSS, i.e., residual sum of squares)
%                         fitMode = 1;
%                         fxn_getDev_LS = @(kernelParams) MC_getDev_CV(paramInd, kernelParams, x_train, yData_train, ifamily, fitMode);
%                         problem_LS = createOptimProblem('fmincon','objective', fxn_getDev_LS,'x0',params0,'lb',lb,'ub',ub,'options',options);
%                         ms_LS = MultiStart('StartPointsToRun', 'bounds', 'Display', 'off');
% %                         params_est_LS = run(ms_LS, problem_LS, nrep);
%                         [params_est_LS] = fxn_parfor(ms_LS, problem_LS, nrep);
%                         %  test
%                         RSS_LS_all(ik, imodel) = MC_getDev_CV(paramInd, params_est_LS, x_test, yData_test, ifamily, fitMode);
%                         
%                     end % end of imodel (parfor)
%                 end % end of ik 
%                 RSS_LS_ave_all(icv, :) = mean(RSS_LS_all,1);
%                 RSS_LS_std_all(icv, :) = std(RSS_LS_all, [], 1);
%                 
%             end % end of icv (cross validation)
%             RSS_LS_ave_allSubj(isubj, :, :) = RSS_LS_ave_all;
%             RSS_LS_std_allSubj(isubj, :, :) = RSS_LS_std_all;
%             [RSS_best, ibest_LS] = min(nanmean(RSS_LS_ave_all));
%             RSS_best_allSubj(isubj) = RSS_best;
%             ibest_LS_allSubj(isubj) = ibest_LS;
%             
%             fprintf('%s done', subjList{isubj})
%             save(sprintf('Data_MC_%s', subjList{isubj}))
%             disp(datetime('now'))
% 
%         end % end of isubj
%         
%         % organize data
%         yData_allSubj_all{icomb, itype} = yData_allSubj;
%         yPred_LS_allSubj_all{icomb, itype} = yPred_LS_allSubj;
%         
%         fprintf('%s vs. %s [%s] done.\n', namesLocComb{combInd(icomb, 1:2)}, namesType{itype})
%         
%     end % end of itype
% end % end of icomb
% 
% fprintf('Model comparison DONE!!\n')
% 
% % save data
% save(MCFileName, '*_allSubj_all', 'xaxis_itp','nmodels', 'namesParams', 'nparams_full', 'ifeature', 'ifamily', ...
%     'paramInd_all', 'ncomb', 'subjList', 'plotAve')
% 
% 
% disp(datetime('now'))
% 
% function [params_est_LS, RSS, flag, output, solution] = fxn_parfor(ms_LS, problem_LS, nrep)
%  [params_est_LS, RSS, flag, output, solution]= run(ms_LS, problem_LS, nrep);
% end