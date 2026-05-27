% 
% % fit the selected model to ALL observers, weighted by number of trials
% clc
% close all
% sampParam0 = @(lb, ub, n) rand(1,n) .* (ub-lb) + lb;
% options = optimoptions('fmincon','MaxIterations',5000,'Display','off');
% nrep = 30; % number of reps to find the global minimum
% nitp = 2; % number of multiples of interpolated pints
% 
% % create functions to calculate AIC/BIC
% % getAIC = @(RSS, nparamsAll, ndata) ndata * log(RSS)+ 2*nparamsAll;
% getAIC_LS = @(RSS, nparamsAll, ndata) ndata * log(RSS) + 2*nparamsAll;
% getAICc_LS = @(RSS, nparamsAll, ndata) ndata * log(RSS) + 2*nparamsAll + 2*nparamsAll*(nparamsAll+1)/(ndata-nparamsAll-1);
% getAIC_ML = @(nLL, nparamsAll, ndata) 2*nLL+ 2*nparamsAll;
% getBIC = @(nLL, nparamsAll, ndata) 2*nLL + nparamsAll * log(ndata);
% 
% % names of params
% namesParams = namesModelParamsAll{ifamily};
% nparams_full = length(namesParams);
% 
% % sizes
% [ncomb, nLocTrue] = size(combInd);
% 
% % x-axis
% if ifeature==2, xaxis_inUse = 2.^xaxis_all{ifeature};
% else, xaxis_inUse = xaxis_all{ifeature};
% end
% nfilters = length(xaxis_inUse);
% xaxis_itp = interp1(1:nfilters, xaxis_inUse, linspace(1, nfilters, nfilters*nitp),'spline');
% 
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
% AIC_LS_allSubj_all = yData_allSubj_all;
% AIC_ML_allSubj_all = yData_allSubj_all;
% AICc_allSubj_all = yData_allSubj_all;
% BIC_allSubj_all = yData_allSubj_all ;
% 
% % estimate time of completion
% ETC = round(12 * ncomb * ntypes * nsubj *nmodels/60/4);
% fprintf('Model comparison starts\n%s [M%d-%s], nmodels=%d. ETC: %d minutes.\n', namesFeature{ifeature}, ifamily, namesFamily_all{ifamily}, nmodels, ETC)
% disp(datetime('now'))
% 
% for icomb = 1%:ncomb
%     
%     % decide the Loc index of the pair to be compared
%     iLocTrue = combInd(icomb, :);
%     
%     for itype = 1:ntypes
%         
%         % extract y data
%         yData_allSubj = squeeze(margK_allSubj{ifeature}(:, itype, icomb, :, :));
%         
%         % empty containers
%         yPred_LS_allSubj = nan(nsubj, nmodels, 2, nfilters*nitp);
%         yPred_ML_allSubj = yPred_LS_allSubj;
%         AIC_LS_allSubj = nan(nsubj, nmodels);
%         AIC_ML_allSubj = nan(nsubj, nmodels);
%         AICc_allSubj = AIC_LS_allSubj;
%         BIC_allSubj = AIC_LS_allSubj;
%         
%         % fit each subj
%         for isubj = 1:nsubj
%             
%             % extract data
%             yData_perSubj = squeeze(yData_allSubj(isubj, :, :));
%             
%             %  empty containers
%             yPred_LS_allModels = nan(nmodels, 2, nfilters*nitp);
%             yPred_ML_allModels = nan(nmodels, 2, nfilters*nitp);
%             AIC_ML_allModels = nan(1,nmodels);
%             AIC_LS_allModels = nan(1,nmodels);
%             AICc_allModels = AIC_ML_allModels;
%             BIC_allModels = AIC_ML_allModels;
%             
%             parfor imodel = 1:nmodels % loop thru param combinations
%                 % the number of data points in one fitting
%                 %                 tic
%                 ndata_itp = length(xaxis_itp) * nLocTrue;
%                 
%                 paramInd = paramInd_all(imodel, :);
%                 nparamsAll = sum(nLocTrue.^paramInd);
%                 
%                 ub = [];
%                 lb = [];
%                 params0 = [];
%                 
%                 % organize ub, lb and param0 given a certain param combination
%                 rng default % For reproducibility
%                 for iparam = 1:nparams_full
%                     if paramInd_all(imodel, iparam) == 1
%                         ub = [ub, ones(1, nLocTrue)*ub_full(iparam)];
%                         lb = [lb, ones(1,nLocTrue)*lb_full(iparam)];
%                         params0 = [params0, sampParam0(lb_full(iparam), ub_full(iparam), nLocTrue)]; % set param0
%                     else
%                         ub = [ub, ub_full(iparam)];
%                         lb = [lb, lb_full(iparam)];
%                         params0 = [params0, sampParam0(lb_full(iparam), ub_full(iparam), 1)]; % set param0
%                     end
%                 end
%                 
%                 assert(length(ub) == nparamsAll)
%                 assert(length(lb) == nparamsAll)
%                 assert(length(params0) == nparamsAll)
%                 
%                 % weight error terms by trial number
%                 %                 if nsubj>1, ntrialsProp = ntrials_allSubj'/sum(ntrials_allSubj); else, ntrialsProp = 1; end
%                 
%                 % fit model ML method (minimizing RSS, i.e., residual sum of squares)
%                 fitMode = 1;
%                 fxn_getDev_LS = @(kernelParams) MC_getDev_IDVD(paramInd, kernelParams, xaxis_itp, yData_perSubj, ifamily, fitMode, nparams_full);
%                 problem_LS = createOptimProblem('fmincon','objective', fxn_getDev_LS,'x0',params0,'lb',lb,'ub',ub,'options',options);
%                 ms_LS = MultiStart('StartPointsToRun', 'bounds', 'Display', 'off');
%                 [params_est_LS, RSS, flag, output, solution] = fxn_parfor(ms_LS, problem_LS, nrep)
%                 
%                 %  organize for plotting
%                 yPred_LS_allModels(imodel, :, :) = RC_predKernel_MC(paramInd, xaxis_itp, yData_perSubj, nparams_full, params_est_LS, ifamily);
%                 AIC_LS_allModels(imodel) = getAIC_LS(RSS, nparamsAll, ndata_itp);
%                 AICc_allModels(imodel) = getAICc_LS(RSS, nparamsAll, ndata_itp);
%                 
%                 % fit model using the ML method (minimizing nLL)
%                 if flagML
%                     fitMode = 2;
%                     fxn_getDev_ML = @(kernelParams) MC_getDev_IDVD(paramInd, kernelParams, xaxis_itp, yData_perSubj, ifamily, fitMode, nparams_full);
%                     problem_ML = createOptimProblem('fmincon','objective', fxn_getDev_ML,'x0',params0,'lb',lb,'ub',ub,'options',options);
%                     ms_ML = MultiStart('StartPointsToRun', 'bounds', 'Display', 'off');
%                     [params_est_ML, nLL, flag, output, solution] = run(ms_ML, problem_ML, nrep);
%                     %  organize for plotting
%                     yPred_ML_allModels(imodel, :, :) = RC_predKernel_MC(paramInd, xaxis_itp, yData_perSubj, nparams_full, params_est_ML, ifamily);
%                     BIC_allModels(imodel) = getBIC(nLL, nparamsAll, ndata_itp);
%                     AIC_ML_allModels(imodel) = getAIC_ML(nLL, nparamsAll, ndata_itp);
%                 end
%                 %             MC_plotSolutions
%                 %                 toc
%             end % end of imodel (parfor)
%             
%             % organize data
%             yPred_LS_allSubj(isubj, :, :, :) = yPred_LS_allModels;
%             AIC_LS_allSubj(isubj, :) = AIC_LS_allModels;
%             AICc_allSubj(isubj, :) = AICc_allModels;
%             if flagML
%                 yPred_ML_allSubj(isubj, :, :, :) = yPred_ML_allModels;
%                 AIC_ML_allSubj(isubj, :) = AIC_ML_allModels;
%                 BIC_allSubj(isubj, :) = BIC_allModels;
%             end
%         end % end of isubj
%         
%         % organize data
%         yData_allSubj_all{icomb, itype} = yData_allSubj;
%         yPred_LS_allSubj_all{icomb, itype} = yPred_LS_allSubj;
%         AIC_LS_allSubj_all{icomb, itype} = AIC_LS_allSubj;
%         AICc_allSubj_all{icomb, itype} = AICc_allSubj;
%         if flagML
%             yPred_ML_allSubj_all{icomb, itype} = yPred_ML_allSubj;
%             AIC_ML_allSubj_all{icomb, itype} = AIC_ML_allSubj;
%             BIC_allSubj_all{icomb, itype} = BIC_allSubj;
%         end
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