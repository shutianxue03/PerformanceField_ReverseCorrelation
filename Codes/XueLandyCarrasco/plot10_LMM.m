
for im = 1%1:2 % 1=CS; 2=pA
    if im == 1, y = cs_allSubj_; else, y = squeeze(metrics_allSubj_(:, :, :, 6));end
    for ifeature = 2%1:2
        fprintf('\n%s: ', namesFeature{ifeature})
        if ifeature == 1
            namesParams = {'PeakAmp', 'Bandwith', 'Baseline'};
            x_med_allSubj = getCI(params_allSubj{ifeature}(:, :, :, itype, 1:3), 1, 2);% peakAmp, bandwidth, baseline
        else
            %             namesParams = {'PeakAmp', 'Bandwith', 'Baseline', 'peakSF', 'Truncation'};
            namesParams = {'PeakAmp', 'Bandwith', 'Baseline', 'peakSF'};
            if flag_plotOctave
                x_med_allSubj = getCI(params_allSubj{ifeature}(:, :, :, itype, [2:4, 1]), 1, 2); % peakAmp, bandwidth, baseline, peakSF, truncation
            else
                x_med_allSubj = getCI(params_allSubj{ifeature}(:, :, :, itype, [2, size(params_allSubj{ifeature}, 5), 4, 1]), 1, 2); % peakAmp, bandwidth, baseline, peakSF, truncation
            end
        end
        nparams = length(namesParams);
        
        % shift vectors with negative values
        
        x_med_allSubjshifted = x_med_allSubj;
        
        for ip = 1:nparams
            x_med_allSubj_ = x_med_allSubj(:, :, ip);
%             if (ifeature==2) && (ip == 3) && flag_plotOctave
%                 x_med_allSubj = getCI(octave(:, :, :, itype), 1, 2); % size(octave)=nsubjx xnB x nLoc x ntypes
%             else
%                 x_med_allSubj = x_med_allSubj;
%             end
%             
%                         % specifics for SF
%             if (ifeature == 2) && (iparam == 3) && ~flag_plotOctave, 
%                 pp_allSubj = params_allSubj{ifeature}(:, :, :, itype, end); 
%             end
            
            if (ifeature==2) && (ip == 3) && flag_plotOctave
                if sum(x_med_allSubj_<=0, 'all')>0
                    fprintf('%s ', namesParams{ip})
                    x_med_allSubjshifted(:, :, ip) = x_med_allSubj_ - min(x_med_allSubj_(:));
                else
                    x_med_allSubjshifted(:, :, ip) = x_med_allSubj_;
                end
            else
                if sum(x_med_allSubjshifted(:, :, ip)<=0, 'all')>0
                    fprintf('%s ', namesParams{ip})
                    x_med_allSubjshifted(:, :, ip) = x_med_allSubj_(:, :, ip) - min(x_med_allSubj_(:, :, ip));
                end
            end
        end % ip
        assert(sum(x_med_allSubjshifted<0, 'all')==0)
        
        % EEI of change in tuningC (x)
        EEI_x = squeeze(getEEI(x_med_allSubjshifted(:, 1, :), x_med_allSubjshifted(:, 2, :)));
        
        % EEI of change in y
        y_med_allSubj = getCI(y, 1, 2);
        EEI_y = getEEI(y_med_allSubj(:, 1), y_med_allSubj(:, 2));
        
        % linear mixed model
        if ifeature==1
            tbl = table(EEI_y, EEI_x(:, 1), EEI_x(:, 2), EEI_x(:, 3), 'VariableNames', {namesM{im}, namesParams{1}, namesParams{2}, namesParams{3}});
            lme = fitlme(tbl, sprintf('%s ~ %s + %s + %s', namesM{im}, namesParams{1}, namesParams{2}, namesParams{3}));
        else
            tbl = table(EEI_y, EEI_x(:, 1), EEI_x(:, 2), EEI_x(:, 3), EEI_x(:, 4), EEI_x(:, 5), 'VariableNames', {namesM{im}, namesParams{1}, namesParams{2}, namesParams{3}, namesParams{4}, namesParams{5}});
            lme = fitlme(tbl, sprintf('%s ~ %s + %s + %s + %s + %s', namesM{im}, namesParams{1}, namesParams{2}, namesParams{3}, namesParams{4}, namesParams{5}));
        end
        EEI_y_pred = fitted(lme);
        
        % beta
        iv=1; lme_text = sprintf('\n%s = %.2f (p = %.3f)\n', namesM{im}, lme.Coefficients{iv,2}, lme.Coefficients{iv,6});
        for iv = 2:nparams+1
            if lme.Coefficients{iv,6} < .001, star = '***';
            elseif lme.Coefficients{iv,6} < .01, star = '**';
            elseif lme.Coefficients{iv,6} < .05, star = '*';
            else, star = '';
            end
            lme_text = [lme_text, sprintf('        + %.2f * %s (p = %.3f)%s\n', lme.Coefficients{iv,2}, lme.Coefficients{iv,1}, lme.Coefficients{iv,6}, star)];
        end
        
        figure('Position', [0 0 500 500]), box on, hold on
        for isubj = 1:nsubj
            plot(EEI_y(isubj), EEI_y_pred(isubj), ['k', markers_allSubj{isubj}])
        end
        axisMin = min([EEI_y_pred; EEI_y]);
        axisMax = max([EEI_y_pred; EEI_y]);
        plot([axisMin, axisMax], [axisMin, axisMax], 'k-')
        xlim([axisMin, axisMax])
        ylim([axisMin, axisMax])
        
        axis square
        xlabel(sprintf('EEI of %s', namesM{im}))
        ylabel(sprintf('Predicted EEI of %s ', namesM{im}))
        title(sprintf('%s - %s vs. %s%s', namesFeature{ifeature}, namesLocComb_{1}, namesLocComb_{2}, lme_text))
        set(findall(gcf, '-property', 'fontsize'), 'fontsize',20)
        
        % save
        folderName = sprintf('%s/%s/LMM/%s/', nameFigFolder, nameEnergySource, nameFileLoc);
        folderDir = dir(folderName); if isempty(folderDir), mkdir(folderName), end
        saveas(gcf, sprintf( '%sn%d_%s_%s.jpg', folderName, nsubj, namesFeature{ifeature}, namesM{im}))
    end % ifeature
end % im
