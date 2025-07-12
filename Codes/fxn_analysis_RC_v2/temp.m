% params_fit_all = cell(1, nparams_full);
% yData = squeeze(margK_allSubj{ifeature}(:, :, itype, iLocTrue, :));
% 
% % Interpolation
% yData_allSubj = nan(nsubj, 2, length(xaxis_itp));
% for isubj=1:nsubj
%     for iLoc = 1:2
%         yData_allSubj(isubj, iLoc, :) = interp1(1:nfilters, squeeze(yData(isubj, iLoc, :)), linspace(1, nfilters, nfilters*nitp),'spline');
%     end
% end
% 
% % empty containers
% %             AIC_allSubj = nan(1, nmodels);
% %             AICc_allSubj = AIC_allSubj;
% %             BIC = AIC_allSubj;
% yPred_SSE_allSubj = nan(nsubj, nmodels, 2, nfilters*nitp);
% AIC_allSubj = nan(nsubj, nmodels);
% iAIC_allSubj = AIC_allSubj;
% AICc_allSubj = AIC_allSubj;
% iAICc_allSubj = AIC_allSubj;
% 
% % fit each subj
% %             figure('Position', [2500 0 1200 200*nsubj])
% for isubj = 1:nsubj
%     
%     %  empty containers
%     yPred_SSE_allModels = nan(nmodels, 2, nfilters*nitp);
%     AIC_allModels = nan(1,nmodels);
%     AICc_allModels = AIC_allModels;
%     
%     for imodel = 1:nmodels % loop thru param combinations
%         % the number of data points in one fitting
%         ndata_itp = length(xaxis_itp) * nLocTrue;
%         
%         %
%         paramInd = paramInd_all(imodel, :);
%         nparamsAll = sum(nLocTrue.^paramInd);
%         
%         ub = [];
%         lb = [];
%         params0 = [];
%         
%         % organize ub, lb and param0 given a certain param combination
%         rng default % For reproducibility
%         for iparam = 1:nparams_full
%             if paramInd_all(imodel, iparam) == 1
%                 ub = [ub, ones(1, nLocTrue)*ub_full(iparam)];
%                 lb = [lb, ones(1,nLocTrue)*lb_full(iparam)];
%                 params0 = [params0, sampParam0(lb_full(iparam), ub_full(iparam), nLocTrue)]; % set param0
%             else
%                 ub = [ub, ub_full(iparam)];
%                 lb = [lb, lb_full(iparam)];
%                 params0 = [params0, sampParam0(lb_full(iparam), ub_full(iparam), 1)]; % set param0
%             end
%         end
%         
%         assert(length(ub) == nparamsAll)
%         assert(length(lb) == nparamsAll)
%         assert(length(params0) == nparamsAll)
%         
%         % weight error terms by trial number
%         %                 if nsubj>1, ntrialsProp = ntrials_allSubj'/sum(ntrials_allSubj); else, ntrialsProp = 1; end
%         
%         % fit model by minimizing SSE
%         fitMode = 1;
%         y_perSubj = squeeze(yData_allSubj(isubj, :, :));
%         fxn_minimizeError_SSE = @(kernelParams) MC_minimizeError_IDVD(paramInd, kernelParams, xaxis_itp, y_perSubj, imodelKernel, fitMode, nparams_full);
%         problem_SSE = createOptimProblem('fmincon','objective', fxn_minimizeError_SSE,'x0',params0,'lb',lb,'ub',ub,'options',options);
%         ms_SSE = MultiStart('StartPointsToRun', 'bounds', 'Display', 'off', 'UseParallel', true);
%         [params_est_SSE, SSE, flag, output, solution] = run(ms_SSE, problem_SSE, nrep);
%         
%         %  organize for plotting
%         yPred_SSE_allModels(imodel, :, :) = RC_predKernel_MC(paramInd, xaxis_itp, y_perSubj, nparams_full, params_est_SSE, imodelKernel);;
%         AIC_allModels(imodel) = getAIC(SSE, nparamsAll, ndata_itp);
%         AICc_allModels(imodel) = getAICc(SSE, nparamsAll, ndata_itp);
%         
%         %             MC_plotSolutions
%         
%     end % end of imodel
%     
%     %% plot selected (5 best) models for each subj
%     %                 buffer = .05;
%     %                 [AIC_sort, iAIC_sort] = sort(AIC_allModels);
%     %                 IC_sort = AIC_sort;
%     %                 iIC_sort = iAIC_sort;
%     %
%     %                 [AICc_sort, iAICc_sort] = sort(AICc_allModels);
%     % %                 ICc_sort = AICc_sort;
%     % %                 iICc_sort = iAICc_sort;
%     %
%     %                 % ========
%     %                 plot_temp
%     %                 % ========
%     
%     % organize data
%     yPred_SSE_allSubj(isubj, :, :, :) = yPred_SSE_allModels;
%     AIC_allSubj(isubj, :) = AIC_allModels;
%     iAIC_allSubj(isubj, :) = iAIC_sort;
%     AICc_allSubj(isubj, :) = AICc_allModels;
%     iAICc_allSubj(isubj, :) = iAICc_sort;
%     
% end % end of isubj
% 
% % organize
% yData_allSubj_all{icomb, itype} = yData_allSubj;
% yPred_SSE_allSubj_all{icomb, itype} = yPred_SSE_allSubj;
% AIC_allSubj_all{icomb, itype} = AIC_allSubj;
% iAIC_s_allSubj_all{icomb, itype} = iAIC_allSubj;
% AICc_allSubj_all{icomb, itype} = AICc_allSubj;
% iAICc_s_allSubj_all{icomb, itype} = iAICc_allSubj;
% 
% fprintf('%s vs. %s [%s] done.\n', namesLocComb{combInd(icomb, 1:2)}, namesType{itype})
