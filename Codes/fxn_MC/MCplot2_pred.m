
ylimits = {[-.05, .16], [-.03, .1]};
% text_x = [-95, .05];
% text1_y = [.08, .08];
% text2_y = [-.07, 0];

for iIC = 1:nIC_
    iBest_unik = unique(iBest_allSubj(:, iIC));
    [freq, iBest_best] = groupcounts(iBest_allSubj(:, iIC));
    nChosen = length(freq);
    
    dev_groupAVE = squeeze(mean(dev_allSubj(:, :, iIC)));
    dev_groupSEM = squeeze(std(dev_allSubj(:, :, iIC)))/sqrt(nsubj);
    [~, irank_groupAVE] = sort(dev_groupAVE); % sort based on averaged dev (for CV) for IC
    iBest_group = irank_groupAVE(1);
    dd_ = dev_groupAVE-min(dev_groupAVE);
    
    if flag_plot
        figure('Position', [0, 0, 2e3, 2e3]), hold on
        for isubj = 1:nsubj
            subjName = subjList{isubj};
            nblocks = nblocks_allSubj(isubj);
            
            nameKernel = sprintf('Data_OOD/%s%d/%s_kernel_%s_%d_%d.mat', subjName, nblocks, subjName, namePatchMode, nORI, nSF);
            
            %         load(nameKernel, 'margORI_perComb', 'margSF_perComb')
            load(nameKernel, 'kernels2D_perComb')
            
            %%  extract y (i.e., kernels)
            [~, nORI, nSF] = size(kernels2D_perComb{1});
            %     kk = nan(nLoc2, nORI, nSF);
            if ifeature==1, marg = nan(nLoc2, nORI); else, marg = nan(nLoc2, nSF); end
            for iiLoc = 1:nLoc2
                %         yData(iiLoc, :) = marg{iLoc_all(iiLoc)}(itype, :);
                kk = squeeze(kernels2D_perComb{iLocComb_all(iiLoc)}(itype, :, :));
                % cut
                if flag_cutMapping
                    kk = kk(cut_ORI, cut_SF);
                end
                % mirror
                if flag_mirrorMapping
                    indMir = (nORI-1)/2;
                    kk_left = kk(1:indMir, :);
                    kk_right = kk(indMir+2:end, :);
                    kk_mid = kk(indMir+1,:);
                    kk_ave = (kk_left + flip(kk_right))/2;
                    kk_mir = [kk_ave; kk_mid; flip(kk_ave)];
                    kk = kk_mir;
                end
                if ifeature==1, marg(iiLoc, :) = mean(kk, 2);
                else, marg(iiLoc, :) = mean(kk, 1);
                end
            end % iiLoc
            %%
            % each is a 8x1 cell, each cell os ntypes x 29
            iBest_idvd = iBest_allSubj(isubj);
            
            subplot(3,5, isubj), hold on
            if iBest_idvd == iBest_group, box on, set(gca,'linewidth',2); else, box off; end
            switch ifeature
                case 1, axis_tuning_ln = axis_tuning{ifeature}; 
                case 2, axis_tuning_ln = 2.^axis_tuning{ifeature}; 
            end
            
            % extra line
            yline(0, '-', 'color', ones(1,3)*.8);
            xline(ifeature-1, '-', 'color', ones(1,3)*.8);
            
            % ylim
            ylim(ylimits{ifeature})
            
            Rsquared_groupBest = []; Rsquared_idvdBest = [];
            for iline = 1:nLoc2
                iLocComb = iLocComb_all(iline);
                
                % raw kernels
                kk = marg(iline, :);
                plot(axis_tuning{ifeature}, kk, 'o', 'color', colors_comb(iLocComb_all(iline), :))
                
                % pred (group best model)
                params_groupBest = params_est_allSubj{isubj}{iBest_group};
                ypred_groupBest = MC_predKernel(paramInd_all(iBest_group, :), axis_tuning_ln, 2, nparams_full, params_groupBest, ifamily);
                plot(axis_tuning{ifeature}, ypred_groupBest(iline, :), '-', 'color', colors_comb(iLocComb_all(iline), :), 'linewidth', 2)
                
                % pred (idvd best model)
                params_idvdBest = params_est_allSubj{isubj}{iBest_idvd};
                ypred_idvdBest = MC_predKernel(paramInd_all(iBest_idvd, :), axis_tuning_ln, 2, nparams_full, params_idvdBest, ifamily);
                plot(axis_tuning{ifeature}, ypred_idvdBest(iline, :), 'k--', 'linewidth', 1)
                
                % Get R squared
                Rsquared_groupBest(iline) = getR2(kk, ypred_groupBest(iline, :), nan);
                Rsquared_idvdBest(iline) = getR2(kk, ypred_idvdBest(iline, :), nan);
                Rsquared_groupBest_allSubj(isubj, iline) = getR2(kk, ypred_idvdBest(iline, :), nan);
            end % iline
            
            %         if isubj==1, legend({'Data', 'Pred (group best)', 'Pred (idvd best)'}, 'Location', 'northeast'), end
            %
            %         title(sprintf('%s\nIDVD best M%d [%s], dev/IC = %.5f\nGroup best: dev/IC = %.5f', ...
            %             subjName, ...
            %             iBest_idvd, num2str(paramInd_all(iBest_idvd, :)), dev_allSubj(isubj, iBest_idvd, iIC), ...
            %             dev_allSubj(isubj, iBest_group, iIC)))
            %
            % xticks
            if ifeature == 2
                xticks(0:.5:2), xticklabels([1, 1.41, 2, 2.83, 4]), xlim([0, 2])
            else, xticks(-90:45:90), xticklabels(-90:45:90), xlim([-90, 90])
            end
            % title and texts
            title(sprintf('%s\nIDVD best M%d [%s] \nIDVD best R^2 = %.0f%% & %.0f%%\nGroup best R^2 = %.0f%% & %.0f%%', ...
                subjName, iBest_idvd, num2str(paramInd_all(iBest_idvd, :)), Rsquared_idvdBest*100, Rsquared_groupBest*100))
            % [text] dev/IC of th group/idvd best model
%             text(text_x(ifeature), text1_y(ifeature), sprintf('IDVD best M%d [%s], dev/IC = %.5f\nGroup best: dev/IC = %.5f', ...
%                 iBest_idvd, num2str(paramInd_all(iBest_idvd, :)), dev_allSubj(isubj, iBest_idvd, iIC), ...
%                 dev_allSubj(isubj, iBest_group, iIC)))
            % [text] R2
%             text(text_x(ifeature), text2_y(ifeature), ...
%                 sprintf('R^2 (group/idvd):\n%s %d%%/%d%%; %s %d%%/%d%%', ...
%                 namesLocComb{iLocComb_all(1)}, round(Rsquared_groupBest(1)*100), round(Rsquared_idvdBest(1)*100), ...
%                 namesLocComb{iLocComb_all(2)}, round(Rsquared_groupBest(2)*100), round(Rsquared_idvdBest(2)*100)))
        end % isubj
        
        if MCmode == 3
            sgtitle(sprintf('n = %d [ORI%d SF%d]\n%s - Family #%d %s\n%s-%s\nBest group model M%d [%s]', ...
                nsubj, nORI, nSF, namesFeature{ifeature}, ifamily, namesFamily_all{ifamily}, namesMCmode{MCmode}, namesIC{iIC}, iBest_group, num2str(paramInd_all(iBest_group, :))))
        else
            sgtitle(sprintf('n = %d [ORI%d SF%d]\n%s - Family #%d %s\n%s\nBest group model M%d [%s] R^2= %.0f%% & %.0f%%', ...
                nsubj, nORI, nSF, namesFeature{ifeature}, ifamily, namesFamily_all{ifamily}, namesMCmode{MCmode}, ...
                iBest_group, num2str(paramInd_all(iBest_group, :)), mean(Rsquared_groupBest_allSubj)*100))
        end
        
        if isempty(dir(sprintf('%spred/', nameFolder_fig))), mkdir(sprintf('%spred/', nameFolder_fig)), end
        if MCmode==3
            saveas(gcf, sprintf('%spred/n%d_M%d_%s.jpg', nameFolder_fig, nsubj, ifamily, namesIC{iIC}))
        else
            saveas(gcf, sprintf('%spred/n%d_M%d_%s.jpg', nameFolder_fig, nsubj, ifamily, namesMCmode{MCmode}))
        end
    end % if flag_plot
end % iIC