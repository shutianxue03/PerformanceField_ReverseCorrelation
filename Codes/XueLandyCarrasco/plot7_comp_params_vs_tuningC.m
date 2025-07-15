
% 2. check the corr among estP/tunC

SX_RC1_setting
close all

labels_estP_ORI = {'Gain', 'Width', 'Baseline'};
labels_estP_SF = {'Peak SF', 'Gain', 'Width', 'Baseline', 'Truncation'};
labels_tunC_ORI = {'PeakAmp', 'Bandwidth', 'Baseline'};
labels_tunC_SF = {'Peak SF', 'PeakAmp', 'Bandwidth', 'Baseline', 'Truncation'};


for ifeature = 1:2
    
    for paramMode = 1:2
        clear margParams
        if paramMode == 1, margParams = estP; xlabels = labels_estP; ylabels = labels_estP; sgtitle_ = 'Estimated params'; figTitle = 'estP';
        else, margParams = tuningC; xlabels = labels_tunC; ylabels = labels_tunC; sgtitle_ = 'Tuning Cs'; figTitle = 'TunC';
        end
        indLoc_corr = repmat(1:nLoc, nsubj, 1);
        indSubj_corr = repmat((1:nsubj)', 1, nLoc);
        
        iplots = 1:nparams*nparams; iplots = reshape(iplots, [nparams, nparams])';
        
        % figure('Position', [0, 0, nparams*300, nparams*300])
        
        for irow = 1:nparams
            % tuning characteristics
            [y, ~, ~, y_neg, y_pos] = getCI(margParams(:, :, :, itype, irow), 1, 2);
            y_ = y-mean(y);
            for icol = 1:nparams
                % estimated params
                [x, ~, ~, x_neg, x_pos] = getCI(margParams(:, :, :, itype, icol), 1, 2);
                x_ = x-mean(x);
                
                % partial corr
                [r_partial, p_partial] = partialcorr([x(:), y(:)], [indLoc_corr(:), indSubj_corr(:)]);
                r_partial = r_partial(2,1); p_partial = p_partial(2,1);
                titles = sprintf('r = %.2f, p = %.3f\n', r_partial, p_partial);
                text_sig = '';
                if p_partial<.05, text_sig = 'sig'; end
                
                % plot
                if irow > icol
                    figure('Position', [0, 0, 500, 500]), hold on
                    nsig = 0;
                    for iLoc = 1:nLoc
                        % sort and connect dots
                        %                 [x_sort, x_i] = sort(x_(:, iLoc));
                        %                 plot(x_sort, y_(x_i, iLoc), '.-', 'color', [colors_comb_(iLoc, :), .4])
                        % idvd dots
                        errorbar(x_(:, iLoc), y_(:, iLoc), x_neg(:, iLoc), x_pos(:, iLoc), 'horizontal', '.','color', colors_comb_(iLoc, :), 'Capsize', 0)
                        errorbar(x_(:, iLoc), y_(:, iLoc), y_neg(:, iLoc), y_pos(:, iLoc), 'vertical', '.','color', colors_comb_(iLoc, :), 'Capsize', 0)
                        
                        for isubj = 1:nsubj
                            plot(x_(isubj, iLoc), y_(isubj, iLoc), markers_allSubj{isubj}, 'markerfacecolor', 'w', 'markeredgecolor', colors_comb_(iLoc, :))
                        end
                        
                        % corr per loc
                        [r_, p_] = corr(x_(:, iLoc), y_(:, iLoc));
                        
                        if p_<.05
                            titles = [titles, sprintf('%s%s', namesLocComb_{iLoc}, getString_starts(p_))];
                            nsig = nsig+1;
                        end
                    end % iLoc
                    
                    % linear regression
                    lm = polyfit(x_(:), y_(:), 1);
                    x_lm2 = linspace(min(x_(:)) - std(x_(:)), max(x_(:)) + std(x_(:)), 2);
                    yfit = polyval(lm, x_lm2);
                    plot(x_lm2, yfit, 'k', 'linewidth', 2, 'handlevisibility', 'off');
                    
                    %
                    xlabel(sprintf('P%d %s', icol, xlabels{icol}))
                    ylabel(sprintf('P%d %s', irow, ylabels{irow}))
                    axis square
                    title(titles)
                    set(findall(gcf, '-property', 'fontsize'), 'fontsize', 18)
                    set(findall(gcf, '-property', 'linewidth'), 'linewidth', 1.2)
                    
                    % save
                    folderName = sprintf('VSS2023/fig/comp_params_tuningC/%s/all%d/', namesFeature{ifeature}, nLoc);
                    folderDir = dir(folderName); if isempty(folderDir), mkdir(folderName), end
                    saveas(gcf, sprintf('%sn%d_%s_p%dp%d_%s%d.jpg', folderName, nsubj, figTitle, icol, irow, text_sig, nsig))
                end % if irow ~= icol
            end % icol
        end% irow
    end % ii
    close all
end % ifeature
