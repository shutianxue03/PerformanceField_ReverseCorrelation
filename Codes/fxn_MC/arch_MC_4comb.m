% 
% % fit the selected model to ALL observers, weighted by number of trials
% 
% sampParam0 = @(lb, ub, n) rand(1,n) .* (ub-lb) + lb;
% options = optimoptions('fmincon','MaxIterations',5000,'Display','off');
% nrep = 30; % number of reps to find the global minimum
% nitp = 2; % number of multiples of interpolated pints
% 
% % create functions to calculate AIC/BIC
% getAIC = @(err, nparamsAll, ndata) ndata * log(err)+ 2*nparamsAll;
% getAICc = @(err, nparamsAll, ndata) ndata * log(err)+ 2*nparamsAll + 2*nparamsAll*(nparamsAll+1)/(ndata-nparamsAll-1);
% getBIC = @(nLL, nparamsAll, ndata) 2*nLL + nparamsAll * log(ndata);
% 
% % names of params
% namesParams = namesModelParamsAll{imodelKernel};
% nparams_full = length(namesParams);
% 
% % file name of the saved estimations
% MCFileName = sprintf('Data/MC_n%d_%s_M%d.mat', nsubj , namesFeature{ifeature}, imodelKernel);
% MCFileDir = dir(MCFileName);
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
% if (imodelKernel == 3) || (imodelKernel == 7) % the truncation term is always left free
%     nmodels = 2^(nparams_full-1);
% else
%     nmodels = 2^nparams_full;
% end
% 
% % empty containers
% AIC_all = cell(ncomb, ntypes);
% AICc_all = AIC_all;
% BIC_all = AIC_all;
% AIC_delta_all = AIC_all;
% AICc_delta_all = AIC_delta_all;
% BIC_delta_all = AIC_delta_all;
% 
% if isempty(MCFileDir)
%     % empty containers
%     yData_allSubj_all = cell(2,ntypes, nmodels);
%     yPred_SSE_allSubj_all = yData_allSubj_all;
%     yPred_nLL_allSubj_all = yData_allSubj_all;
%     AIC_allSubj_all =yData_allSubj_all ;
%     AICc_allSubj_all = yData_allSubj_all;
%     BIC_allSubj_all = yData_allSubj_all ;
%     SSE_all = nan(ncomb, ntypes, nmodels);
%     
%     % estimate time of completion
%     ETC = round(.11 * nsubj * ntypes * ncomb * nmodels * nrep/60);
%     fprintf('Model comparison starts: %s [M%d], nmodels=%d. ETC: %d minutes.\n', namesFeature{ifeature}, imodelKernel, nmodels, ETC)
%     disp(datetime('now'))
%     
%     % extract ub and lb
%     ub_full = ub_full_all{imodelKernel};
%     lb_full = lb_full_all{imodelKernel};
%     
%     for icomb = 1:ncomb
%         
%         % decide the index of the comparison pair
%         iLocTrue = combInd(icomb, :);
%         
%         % decide the param index (1 = specified, 0 = shared)
%         if (imodelKernel == 3) || (imodelKernel == 7) % for model with a truncation term, it's always left free to vary
%             paramInd_all_ = fxn_getParamInd(nparams_full-1);
%             paramInd_all = [paramInd_all_, ones(nmodels,1)];
%         else
%             paramInd_all = fxn_getParamInd(nparams_full);
%         end
%         
%         for itype = 1:ntypes
%             params_fit_all = cell(1, nparams_full);
%             yData = squeeze(margK_allSubj{ifeature}(:, :, itype, iLocTrue, :));
%             
%             % Interpolation
%             yData_allSubj = nan(nsubj, 2, length(xaxis_itp));
%             for isubj=1:nsubj
%                 for iLoc = 1:2
%                     yData_allSubj(isubj, iLoc, :) = interp1(1:nfilters, squeeze(yData(isubj, iLoc, :)), linspace(1, nfilters, nfilters*nitp),'spline');
%                 end
%             end
%             
%             % the number of data points in one fitting
%             ndata_itp = length(xaxis_itp) * nLocTrue * nsubj;
%             
%             % empty containers
%             AIC = nan(1, nmodels);
%             AICc = AIC;
%             BIC = AIC;
%             
%             for imodel = 1:nmodels % loop thru param combinations
%                 
%                 %
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
%                 fitMode = 1;
%                 fxn_minimizeError_SSE = @(kernelParams) MC_minimizeError(paramInd, kernelParams, xaxis_itp, yData_allSubj, imodelKernel, fitMode, nparams_full, ntrials_allSubj);
%                 problem_SSE = createOptimProblem('fmincon','objective', fxn_minimizeError_SSE,'x0',params0,'lb',lb,'ub',ub,'options',options);
%                 ms_SSE = MultiStart('StartPointsToRun', 'bounds', 'Display', 'off', 'UseParallel', true);
%                 [params_est_SSE, SSE, flag, output, solution] = run(ms_SSE, problem_SSE, nrep);
%                 SSE_all(icomb, itype, imodel) = SSE;
%                 %             MC_plotSolutions
%                 
%                 % calculate the fitting line & error for each subj
%                 yPred_SSE_allSubj = nan(nsubj, 2, nfilters*nitp);
%                 AICc_allSubj = nan(nsubj,1);
%                 AIC_allSubj = nan(nsubj,1);
%                 for isubj = 1:nsubj
%                     dd = yData_allSubj(isubj, :, :);
%                     yPred_ = RC_predKernel_MC(paramInd, xaxis_itp, squeeze(dd), nparams_full, params_est_SSE, imodelKernel);
%                     yPred_SSE_allSubj(isubj, :, :) = yPred_;
%                     AIC_allSubj(isubj) = getAIC(sumsqr(yPred_(:) - dd(:))/length(dd(:)), nparamsAll, length(dd(:)));
%                     AICc_allSubj(isubj) = getAICc(sumsqr(yPred_(:) - dd(:))/length(dd(:)), nparamsAll, length(dd(:)));
%                 end
%                 
%                 %% use MLE to fit data
% %                 fitMode = 2;
% %                 fxn_minimizeError_nLL = @(kernelParams) MC_minimizeError(paramInd, kernelParams, xaxis_itp, yData_allSubj, imodelKernel, fitMode, nparams_full, ntrials_allSubj);
% %                 problem_nLL = createOptimProblem('fmincon','objective', MC_minimizeError_nLL,'x0',params0,'lb',lb,'ub',ub,'options',options);
% %                 ms_nLL = MultiStart('StartPointsToRun', 'bounds', 'Display', 'off', 'UseParallel', true);
% %                 [params_est_nLL, nLL] = run(ms_nLL, problem_nLL, nrep);
%                 
%                 % calculate the fitting line & error for each subj
%                 yPred_nLL_allSubj = nan(nsubj, 2, nfilters*nitp);
%                 BIC_allSubj = nan(nsubj,1);
% %                 for isubj = 1:nsubj
% %                     dd = yData_allSubj(isubj, :, :);
% %                     yPred_ = RC_predKernel_MC(paramInd, xaxis_itp, squeeze(dd), nparams_full, params_est_nLL, imodelKernel);
% %                     yPred_nLL_allSubj(isubj, :, :) = yPred_;
% %                     BIC_allSubj(isubj) = getBIC(sumsqr(yPred_(:) - dd(:))/length(dd(:)), nparamsAll, length(dd(:)));
% %                 end
%                 
%                 % AIC and BIC
%                 AIC(imodel) = getAIC(SSE, nparamsAll, ndata_itp);
%                 AICc(imodel) = getAICc(SSE, nparamsAll, ndata_itp);
% %                 BIC(imodel) = getBIC(nLL, nparamsAll, ndata_itp);
%                 
%                 % organize data
%                 yData_allSubj_all{icomb, itype, imodel} = yData_allSubj;
%                 yPred_SSE_allSubj_all{icomb, itype, imodel} = yPred_SSE_allSubj;
%                 yPred_nLL_allSubj_all{icomb, itype, imodel} = yPred_nLL_allSubj;
%                 AIC_allSubj_all{icomb, itype, imodel} = AIC_allSubj;
%                 AICc_allSubj_all{icomb, itype, imodel} = AICc_allSubj;
%                 BIC_allSubj_all{icomb, itype, imodel} = BIC_allSubj;
%                 
%             end % end of imodel
%             
%             fprintf('%s vs. %s [%s] done.\n', namesLocComb{combInd(icomb, 1:2)}, namesType{itype})
%             
%             % raw AIC
%             AIC_all{icomb, itype} = AIC;
%             AICc_all{icomb, itype} = AICc;
%             BIC_all{icomb, itype} = BIC;
%             
%             % delta AIC
%             AIC_delta_all{icomb, itype} = AIC - min(AIC);
%             AICc_delta_all{icomb, itype} = AICc - min(AICc);
%             BIC_delta_all{icomb, itype} = BIC - min(BIC);
%             
%         end % end of itype
%     end % end of icomb
%     
%     % save
%     save(MCFileName, '*_all', 'xaxis_itp','nmodels', 'namesParams', 'nparams_full', 'ifeature', 'imodelKernel',...
%         'ncomb', 'subjList', 'plotAve')
% else
%     fprintf('MC file loaded.\n')
%     load(MCFileName)
% end
% 
% %% plot fitting line for each subj
% nIC = 2;
% IC_all_all = {AIC_all, AICc_all, BIC_all};
% IC_type_all = {'AIC', 'AICc', 'BIC'};
% IC_delta_all_all = {AIC_delta_all, AICc_delta_all, BIC_delta_all};
% IC_allSubj_all_all = {AIC_allSubj_all, AICc_allSubj_all, BIC_allSubj_all};
% 
% publishOptionsMC = publishOptions;
% publishOptionsMC.outputDir = 'publishedPDFs/MC/';
% 
% imodel_opt_all = cell(1,3);
% imodel_opt_IC_all = cell(1,3);
% imodel_opt_paramInd_all = cell(1,3);
% 
% for iIC = 1:nIC % AIC,AICc, BIC
%     if iIC == 3
%         yPred_allSubj_all = yPred_nLL_allSubj_all;
%     else
%         yPred_allSubj_all = yPred_SSE_allSubj_all;
%     end
%     
%     IC_all = IC_all_all{iIC};
%     IC_type = IC_type_all{iIC};
%     IC_delta_all = IC_delta_all_all{iIC};
%     IC_allSubj_all = IC_allSubj_all_all{iIC};
%     
%     if plotFlag
%         publish('MCplot_4comb', publishOptionsMC);
%         movefile([publishOptionsMC.outputDir,'MCplot_4comb.pdf'], sprintf('%sMC_n%d_%s_M%d_%s.pdf',  publishOptionsMC.outputDir, nsubj, namesFeature{ifeature}, imodelKernel, IC_type));
%         close all
%     else
%         MC_plotIDVD
%     end
%     
%     imodel_opt_all{iIC} = imodel_opt;
%     imodel_opt_IC_all{iIC} = imodel_opt_IC;
%     imodel_opt_paramInd_all{iIC} = imodel_opt_paramInd;
%     
%     fprintf('DONE\n')
% end
% 
% save(MCFileName, 'imodel_opt_all', 'imodel_opt_IC_all', 'imodel_opt_paramInd_all', '-append')
% 
% %% plot ICs
% markers_IC = {'o', '*', '+'};
% figure('Position', [0, 200, 1000 800])
% 
% for icomb = 1:ncomb
%     
%     for itype = 1:ntypes
%         
%         subplot(ntypes, ncomb, ncomb*(itype-1)+icomb), hold on
%         xticklabels_IC = cell(1,nIC);
%         
%         for iIC=1:nIC
%             
%             imodel_opt_IC = imodel_opt_IC_all{iIC};
%             plot(iIC, imodel_opt_IC(icomb, itype), ['k', markers_IC{iIC}])
%             imodel_opt = imodel_opt_all{iIC};
%             imodel_opt_paramInd = imodel_opt_paramInd_all{iIC};
%             xticklabels__ = sprintf('M%d', imodel_opt(icomb, itype));
%             for ip = 1:nparams_full
%                 xticklabels__ = sprintf('%s\\newline%d', xticklabels__ , imodel_opt_paramInd{icomb, itype}(ip));
%             end
%             xticklabels_IC{iIC} = xticklabels__;
%         end
%         
%         xticks(1:nIC)
%         xticklabels(xticklabels_IC)
%         xlim([.5, nIC+.5])
%         
%         if icomb == 1, ylabel(namesType{itype}), end
%         if itype == 1, title(sprintf('%s vs. %s', namesLocComb{combInd(icomb, [1,2])})), end
%         if icomb + itype == 2, legend(IC_type_all, 'Location', 'best'), end
%     end
% end
% 
% set(findall(gcf, '-property', 'FontSize'), 'FontSize',12)
% sgtitle(sprintf('%s-M%d', namesFeature{ifeature}, imodelKernel), 'fontsize', 25)
% 
% saveas(gcf, sprintf('%s/n%d_%s_M%d.jpg', publishOptionsMC.outputDir, nsubj, namesFeature{ifeature}, imodelKernel))
% 
