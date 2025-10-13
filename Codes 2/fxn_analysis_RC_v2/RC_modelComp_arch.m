% 
% 
% sampParam0 = @(lb, ub, n) rand(1,n) .* (ub-lb) + lb;
% % options = optimoptions('fmincon','MaxIterations',5000,'Display','off', 'PlotFcn', 'optimplotx');
% % options = optimoptions('fmincon','MaxIterations',5000,'Display','off', 'PlotFcn', 'optimplotstepsize');
% options = optimoptions('fmincon','MaxIterations',5000,'Display','off');
% nrep = 30; % number of reps to find the global minimum
% 
% 
% %%
% % frequently used imodelKernel_perF:
% % 1 = gaussian
% % 2 = log parabola, AF
% % 3 =  truncated log parabola
% % 8 = DoG [ORI]
% 
% imodelKernel_perF = [8,3];
% 
% %%
% [n_iLocTrue, nLocTrue] = size(combInd);
% 
% AIC_delta_perF = cell(n_iLocTrue, 2, ntypes);
% AICc_delta_perF = AIC_delta_perF;
% BIC_delta_perF = AIC_delta_perF;
% 
% for ifeature = 1
%     
%     if ifeature==2, xaxis_inUse = 2.^xaxis_all{ifeature}; else, xaxis_inUse = xaxis_all{ifeature}; end
%     nfilters = length(xaxis_inUse);
%     imodelKernel = imodelKernel_perF(ifeature);
%     
%     namesParams = namesModelParamsAll{imodelKernel};
%     nparams_full = length(namesParams);
%     if imodelKernel == 3 % the truncation term is always left free
%         nmodels = 2^(nparams_full-1);
%     else
%         nmodels = 2^nparams_full;
%     end
%     
%     % estimate time taken
%     fprintf('Model comparison starts. ETA: %d minutes.\n', round(.16 * ntypes * 4 * nmodels * nrep/60))
%     disp(datetime('now'))
%     
%     % decide the grid of the plots
%     nrow = ceil(sqrt(nmodels));
%     ncol = nmodels/nrow;
%     assert(nrow*ncol >= nparams_full)
%     
%     % extract ub and lb
%     ub_full = ub_full_all{imodelKernel};
%     lb_full = lb_full_all{imodelKernel};
%     
%     for icomb = 1:n_iLocTrue
%         % decide the index of the comparison pair
%         iLocTrue = combInd(icomb, :);
%         
%         % for model with a truncation term, it's always left free to vary
%         if imodelKernel == 3
%             paramInd_all = fxn_getParamInd(nparams_full-1);
%             paramInd_all = [paramInd_all, ones(nmodels,1)]; 
%         else
%             paramInd_all = fxn_getParamInd(nparams_full);
%         end
%         
%         
%         % functions to calculate AIC/BIC
%         getAIC = @(err, nparamsAll, ndata) ndata * log(err)+ 2*nparamsAll;
%         getAICc = @(err, nparamsAll, ndata) ndata * log(err)+ 2*nparamsAll + 2*nparamsAll*(nparamsAll+1)/(ndata-nparamsAll-1);
%         getBIC = @(nLL, nparamsAll, ndata) 2*nLL + nparamsAll * log(ndata);
%         
%         for itype = 1:ntypes
% 
%             params_fit_all = cell(1, nparams_full);
%             yData = squeeze(margK_perComb{ifeature}(itype, iLocTrue, :));
%             
%             % interpolate
% %             nitp = 50;
% %             xaxis_itp = interp1(1:nfilters, xaxis_inUse, 1:.5:nfilters,'spline');
% %             yData_itp = nan(2,length(xaxis_itp));
% %             yData_itp(1,:) = interp1(1:nfilters, yData(1,:), 1:.5:nfilters,'spline');
% %             yData_itp(2,:) = interp1(1:nfilters, yData(2,:), 1:.5:nfilters,'spline');
%             
%             % not 
%             xaxis_itp = xaxis_inUse;
%             yData_itp = yData;
%             
%             % the number of data points
%             ndata_itp = length(xaxis_itp) * nLocTrue;
%             
%             AIC = nan(1, nmodels);
%             AICc = nan(1, nmodels);
%             BIC = nan(1, nmodels);
%             
%             if plotFlag, figure('Position', [0 400 ncol*200 nrow*200]), end
%             
%             for imodel = 1:nmodels % loop thru param combinations
%                 
%                 modelInd = paramInd_all(imodel, :);
%                 nparamsAll = sum(nLocTrue.^modelInd);
%                 ub = [];
%                 lb = [];
%                 params0 = [];
%                 
%                 % organize ub, lb and param given a certain param combination
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
%                 %%%%%%%%%%%%%%
%                 %    use SSE to fit data    %
%                 %%%%%%%%%%%%%%
%                 fitMode = 1;
%                 kernelFit_fxn_SSE = @(kernelParams) kernelFit_modelComp(modelInd, kernelParams, xaxis_itp, yData_itp, imodelKernel, fitMode, nparams_full);
%                 problem = createOptimProblem('fmincon','objective',kernelFit_fxn_SSE,'x0',params0,'lb',lb,'ub',ub,'options',options);
%                 ms = MultiStart('StartPointsToRun', 'bounds', 'Display', 'off');
%                 [kernelParams_est_SSE, SSE] = run(ms, problem, nrep);
% %                 [kernelParams_est_SSE, SSE, flag, outpt, mins_all] = run(ms, problem, nrep);
% 
%                 y_pred_SSE = RC_predKernel_modelComp(modelInd, xaxis_itp, yData_itp, nparams_full, kernelParams_est_SSE, imodelKernel);
%                 
%                 %%%%%%%%%%%%%%
%                 %    use MLE to fit data    %
%                 %%%%%%%%%%%%%%
% %                 fitMode = 2;
% %                 kernelFit_fxn_nLL = @(kernelParams) kernelFit_modelComp(modelInd, kernelParams, xaxis_itp, yData_itp, imodelKernel, fitMode, nparams_full);
% %                 [kernelParams_est_nLL, nLL, exitflag, output, lambda,grad,hessian] = fmincon(kernelFit_fxn_nLL, params0, [], [], [], [], lb, ub, [], options);
% %                 y_pred_nLL = RC_predKernel_modelComp(modelInd, xaxis_itp, yData_itp, nparams_full, kernelParams_est_nLL, imodelKernel);
%                 
%                 %%%%%%%%%%%%
%                 %       AIC and BIC      %
%                 %%%%%%%%%%%%
% %                 AIC(imodel) = getAIC(SSE, nparamsAll, ndata_itp);
%                 AICc(imodel) = getAICc(SSE, nparamsAll, ndata_itp);
% %                 BIC(imodel) = getBIC(nLL, nparamsAll, ndata_itp);
%                 
%                 %%%%%%%%
%                 %       plot       %
%                 %%%%%%%%
%                 printFlag = 0; % 0=do NOT print the estimated params
%                 xaxis = xaxis_itp;
%                 yData = yData_itp;
%                 y_pred = y_pred_SSE; % could change to y_pred_nLL if needed
%                 kernelParams_est = kernelParams_est_SSE; % could change to kernelParams_est_nLL if needed
%                 if plotFlag, fxn_plotTuningComb_MC, end
%                 
%             end
%             
% %             fprintf('%s vs. %s [%s] done.\n', namesLocComb{combInd(icomb, 1:2)}, namesType{itype})
% 
%             AIC_delta = AIC - min(AIC);
%             AICc_delta = AICc - min(AICc);
%             BIC_delta = BIC - min(BIC);
%             
%             AIC_delta_perF{icomb, ifeature, itype} = AIC_delta;
%             AICc_delta_perF{icomb, ifeature, itype} = AICc_delta;
%             BIC_delta_perF{icomb, ifeature, itype} = BIC_delta;
%             
%             % plot
%             if plotFlag
%                 title_ = 'AIC_c'; 
%                 delta_allSubj2 = AICc_delta;
%                 subplot(nrow+1, ncol, ncol*nrow+1: ncol*(nrow+1)), hold on
%                 plotIDVD = 0; 
%                 flagLongXTicks = 1;% show the param index in x-axis(1 for free and 0 for fixed)
%                 % ==================
%                 RCplot_AICc
%                 % ==================
%                 ylabel(title_)
%                 set(findall(gcf, '-property', 'FontSize'), 'FontSize',12)
%                 
%                 % highlight the optimal model
%                 [~, imodel_opt] = min(AICc_delta);
%                 set(ax_all(imodel_opt),'color',[1 .8 .75])
%                 set(ax_all(imodel_opt),'box','on')
%             end
%             
%             % plot all the fitted params
% %             plotFittedParams_MC
% 
%         end
%     end
% end
% 
% 
% 
% 
