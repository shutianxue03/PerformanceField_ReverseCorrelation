
% params barplot
load params_n4_loc5 % params_allLoc_loc5: 5xnparams
load params_n4_PF3 % params_allLoc_PF3: 3xnparams, HM, VM, peri

idata_s_all = [1,8;6,7;5,3];
ncomp = size(idata_s_all, 1);

itype = 3;

ylimits = cell(2,3);

% ylimits of SF
nsteps = 3;
yticks_peak = [1.5, 1.8, 2.1];
yticks_gain_SF = [0, .06, .12];
yticks_width_SF = [0, .4 .8];

ylimits{1,1} = [yticks_peak; yticks_gain_SF; yticks_width_SF;  linspace(-.005, 0, length(yticks_gain_SF))];     % center vs. periphery
ylimits{1,2} = [yticks_peak; yticks_gain_SF; yticks_width_SF; linspace(-.005, .005, length(yticks_gain_SF))]; % HM vs. VM
ylimits{1,3} = [yticks_peak; yticks_gain_SF; yticks_width_SF; linspace(-.01, 0, length(yticks_gain_SF))];        % LVM vs. UVM

% ylimits of ori
yticks_gain_ori = 0:3:6;
yticks_width_ori = 0:10:20;
ylimits{2,1} = [yticks_gain_ori; yticks_width_ori;  linspace(-.005, 0, length(yticks_gain_ori))];     % center vs. periphery
ylimits{2,2} = [yticks_gain_ori; yticks_width_ori; linspace(-.005, .005, length(yticks_gain_ori))]; % HM vs. VM
ylimits{2,3} = [yticks_gain_ori; yticks_width_ori; linspace(-.01, 0, length(yticks_gain_ori))];        % LVM vs. UVM

%% plot the kernel
% params_allComp_both = cell(1,2);
%
% for icomp = 1:ncomp
%
%     % extract kernel for comparison
%     kernel2D_comp1 = kernels2D_allSubj_comp{idata_s_all(icomp, 1)};
%     kernel2D_comp2 = kernels2D_allSubj_comp{idata_s_all(icomp, 2)};
%
%     for ikernel = 1:2%:2
%
%         switch ikernel
%             case 1
%                 nparams = nparamsSF;
%                 paramsNames = modelParamsNamesSF;
%                 modelInd = modelSF;
%                 idim = 3;
%                 nfilters = nfiltersSF_cut;
%                 xaxis = filtersSF_all_log_cut;
%                 params_allSubj = paramsSF_allSubj;
%             case 2
%                 nparams = nparamsOri;
%                 paramsNames = modelParamsNamesOri;
%                 modelInd = modelOri;
%                 idim = 4;
%                 nfilters = nfiltersOri;
%                 xaxis = filtersOri_all-90;
%                 params_allSubj = paramsOri_allSubj;
%         end
%
%         params_allSubj_allComp = nan(nsubj, nparams);
%
%         % marginalize for all subjects
%         kernels_m_allSubj1 = squeeze(mean(kernel2D_comp1, idim));
%         kernels_m_allSubj2 = squeeze(mean(kernel2D_comp2, idim));
%
%         % get the mean and errorband of tuning fxn
%         if nsubj > 1
%             % get weight-average based on ntrials
%             kernels_m_ave1 = 0;
%             kernels_m_ave2 = 0;
%
%             for isubj = 1:nsubj
%                 kernels_m_ave1 = kernels_m_ave1 + squeeze(kernels_m_allSubj1(isubj, itype, :)) * ntrialsProp(isubj);
%                 kernels_m_ave2 = kernels_m_ave2 + squeeze(kernels_m_allSubj2(isubj, itype, :)) * ntrialsProp(isubj);
%             end
%
%             % get error
%             kernels_err1 = squeeze(nanstd(kernels_m_allSubj1(:,itype, :), [], 1)) / sqrt(nsubj);
%             kernels_err2 = squeeze(nanstd(kernels_m_allSubj2(:,itype, :), [], 1)) / sqrt(nsubj);
%         else
%             kernels_m_ave1 = squeeze(kernels_m_allSubj1);
%             kernels_m_ave2 = squeeze(kernels_m_allSubj2);
%             kernels_err1 = zeros(1, nfilters);
%             kernels_err2 = kernels_err1;
%         end
%
%         % fit for each subj (to get errorbar of params)
%         for isubj = 1:nsubj
%             % extract data for each observer and itype
%             kernels_m_1 = squeeze(kernels_m_allSubj1(isubj, itype, :));
%             % fit, get ypred and params
%             if ikernel==1, [kernels_pred, params, Rsquared, ~] = SX_sim08_fit(2.^xaxis, reshape(kernels_ave_, 1, length(xaxis)), modelSF, fitMode);
%             else, [kernels_pred, params, Rsquared, ~] = SX_sim08_fit(xaxis, reshape(kernels_ave_, 1, length(xaxis)), modelOri, fitMode);
%             end
%             params_allSubj_allComp(isubj, :) = params;
%         end
%
%         % plot
%         figure, hold on
%
%
%
%
%     end
% end


%% plot the params
for ikernel = 1%:2
    switch ikernel
        case 1
            nparams = nparamsSF;
            paramsNames = modelParamsNamesSF;
            modelInd = modelSF;
            idim = 3;
            nfilters = nfiltersSF_cut;
            xaxis = filtersSF_all_log_cut;
            params_allSubj = paramsSF_allSubj;
        case 2
            nparams = nparamsOri;
            paramsNames = modelParamsNamesOri;
            modelInd = modelOri;
            idim = 4;
            nfilters = nfiltersOri;
            xaxis = filtersOri_all-90;
            params_allSubj = paramsOri_allSubj;
    end
    
    nparams_cut = nparams - 2; % do not plot baseline
    
    if publishFlag, figure('Position', [0 200 ncomp*300 (nparams_cut)*200]), end
    for icomp = 1:ncomp
        if ~publishFlag, figure('Position', [0 200 300 (nparams_cut)*200]), end
        
        idata_s = idata_s_all(icomp, :);
        
        for iparam = 1:nparams_cut
            if publishFlag, subplot(nparams_cut, ncomp, (iparam-1)*ncomp+icomp), hold on, box on, ax = gca; ax.YGrid = 'on';
            else, subplot(nparams_cut, 1, iparam), hold on, box on, ax = gca; ax.YGrid = 'on';
            end
            
            params_bysubj_forcomp = cell(1,2);
            for idata = 1:2 % the two data set to be compared
                
                % estimate params of each subj, then average
                params_bysubj = nan(nsubj, nparams);
                kernels2D_s = kernels2D_allSubj_comp{idata_s(idata)};
                if nsubj>1
                    for isubj = 1:nsubj
                        k_m = mean(kernels2D_s(isubj, itype, :, :), idim);
                        if ikernel == 1, [kernels_pred, params, Rsquared, ~] = SX_sim08_fit(2.^xaxis, reshape(k_m, 1, nfilters), modelInd, fitMode);
                        else, [kernels_pred, params, Rsquared, ~] = SX_sim08_fit(xaxis, reshape(k_m, 1, nfilters), modelInd, fitMode);
                        end
                        params_bysubj(isubj, :) = params;
                    end
                else
                    k_m = squeeze(mean(kernels2D_s(itype, :,:),idim-1));
                    if ikernel == 1, [kernels_pred, params, Rsquared, ~] = SX_sim08_fit(2.^xaxis, reshape(k_m, 1, nfilters), modelInd, fitMode);
                    else, [kernels_pred, params, Rsquared, ~] = SX_sim08_fit(xaxis, reshape(k_m, 1, nfilters), modelInd, fitMode);
                    end
                    params_bysubj = params;
                end
                params_bysubj_forcomp{idata} = params_bysubj;
                
                % plot bar
                bar(idata, mean(params_bysubj(:, iparam)), 'barwidth', .3, 'facecolor', colors_comp(idata_s(idata), :), 'edgecolor', colors_comp(idata_s(idata), :))
                if (ikernel == 1) && (idata ==1), yline(2); end
                % plot error
                if nsubj>1, errorbar(idata, mean(params_bysubj(:, iparam)), std(params_bysubj(:, iparam))/sqrt(nsubj), 'k'), end
            end
            
            xticks([1,2])
            xticklabels(locCompNames(idata_s))
            ax = gca; ax.FontSize = sz_ticks;
            xlim([.5, 2.5])
            yticks(round(ylimits{ikernel, icomp}(iparam, :), 2))
            ylim(ylimits{ikernel, icomp}(iparam, [1,end]))
            
            if iparam == 1, ylabel(sprintf('%s (a.u.)', paramsNames{iparam}))
            else, ylabel(sprintf('%s (deg)', paramsNames{iparam}))
            end
            title(paramsNames{iparam}, 'fontsize', sz_title)
            
            % paired ttest
            fprintf('\n%s vs. %s:\n', locCompNames{idata_s(1)}, locCompNames{idata_s(2)})
            if nsubj>1
                a = params_bysubj_forcomp{1}(:, iparam);
                b = params_bysubj_forcomp{2}(:, iparam);
                [~,p,CI,stats] = ttest(a,b);
                fprintf('%s: t(%d) = %.3f, p = %.3f\n', paramsNames{iparam}, stats.df, stats.tstat, p)
            end
        end
        
        set(findall(gcf, '-property', 'linewidth'), 'linewidth',2)
        
        % save figure
        if ~publishFlag
            if nsubj>1, saveas(gcf, sprintf('publishedPDFs/fig/n%d_params_%s_%s_%s.jpg', nsubj, kernelNames{ikernel}, locCompNames{idata_s(1)}, locCompNames{idata_s(2)}))
            else, saveas(gcf, sprintf('publishedPDFs/fig/%s_params_%s_%s_%s.jpg', subjName, kernelNames{ikernel},locCompNames{idata_s(1)}, locCompNames{idata_s(2)}))
            end
        end
    end
end
