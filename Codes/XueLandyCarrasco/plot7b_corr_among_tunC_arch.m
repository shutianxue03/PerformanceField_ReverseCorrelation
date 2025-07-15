%  check the corr among estP/tunC

% SX_RC1_setting
clc
close all

labels_estP_ORI = namesTunC_unit_perF{ifamily_perF(1), 1};
labels_estP_SF = namesTunC_unit_perF{ifamily_perF(2), 1};
labels_tunC_ORI = namesTunC_unit_perF{ifamily_perF(1), 2};
labels_tunC_SF = namesTunC_unit_perF{ifamily_perF(2), 2};

sz_marker=20;

% if ifeature==1, nparams = 3; else, nparams = 4; end

for paramMode = 2
    clear margParams
    %     if paramMode == 1
    %         figTitle = 'estP';
    %         if ifeature==1, p = estP_ORI_allSubj; labels = labels_estP_ORI;% ylabels = labels_estP;
    %         else, p = estP_SF_allSubj; labels = labels_estP_SF; %ylabels = labels_estP;
    %         end
    %     else, figTitle = 'TunC';
    %         if ifeature==1, p = tuningC_ORI_allSubj; labels = labels_tunC_ORI;% ylabels = labels_estP;
    %         else, p = tuningC_SF_allSubj; labels = labels_tunC_SF; %ylabels = labels_estP;
    %         end
    if ifeature==1, tunC = tunC_ORI; labels = labels_tunC_ORI;% ylabels = labels_estP;
    else, tunC = tunC_SF; labels = labels_tunC_SF; %ylabels = labels_estP;
    end
    
    %     end
    nparams = length(labels);
    nLoc = size(tunC, 3);
    indLoc_corr = repmat(1:nLoc, nsubj, 1);
    indSubj_corr = repmat((1:nsubj)', 1, nLoc);
    
    iplots = 1:nparams*nparams; iplots = reshape(iplots, [nparams, nparams])';
    
    for irow = 1:nparams
        % param_y
        [y, ~, ~, y_neg, y_pos] = getCI(tunC(:, :, :, itype, irow), 1, 2);
        y_ = y-mean(y);
        for icol = 1:nparams
            % param_x
            [x, ~, ~, x_neg, x_pos] = getCI(tunC(:, :, :, itype, icol), 1, 2);
            x_ = x-mean(x);
            
            % partial corr
            %             [r_partial, p_partial] = partialcorr([x(:), y(:)], [indLoc_corr(:), indSubj_corr(:)]);
            [r_partial, p_partial] = partialcorr([x(:), y(:)], indLoc_corr(:));
            r_partial = r_partial(2,1); p_partial = p_partial(2,1);
            titles = sprintf('r = %.2f, p = %.3f\n', r_partial, p_partial);
            text_sig = '';
            if p_partial<.05, text_sig = 'sig'; end
            
            % plot
            if irow > icol % so that the current pair is not repetitive
                fprintf('[%d, %d] %s vs. %s\n', irow, icol, labels{irow}, labels{icol})
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
                        plot(x_(isubj, iLoc), y_(isubj, iLoc), markers_allSubj{isubj}, ...
                            'markerfacecolor', 'w', 'markeredgecolor', colors_comb_(iLoc, :), 'markersize', sz_marker)
                    end
                    
                    % corr per loc
                    [r_, p_] = corr(x_(:, iLoc), y_(:, iLoc));
                    titles = [titles, sprintf('%s: r=%.2f, p=%.2f\n', namesLocComb_{iLoc}, r_, p_)];
                    if p_<.05, nsig = nsig+1; end
                end % iLoc
                
                % linear regression
                lm = polyfit(x_(:), y_(:), 1);
                x_lm2 = linspace(min(x_(:)) - std(x_(:)), max(x_(:)) + std(x_(:)), 2);
                yfit = polyval(lm, x_lm2);
                plot(x_lm2, yfit, 'k', 'linewidth', 2, 'handlevisibility', 'off');
                
                %
                axis square
                %                     xlabel(sprintf('P%d %s', icol, labels{icol}))
                %                     ylabel(sprintf('P%d %s', irow, labels{irow}))
                xlabel(sprintf('%s %s', namesFeature{ifeature}, labels{icol}))
                ylabel(sprintf('%s %s', namesFeature{ifeature}, labels{irow}))
                title(sprintf('[%s] %s', namesFeature{ifeature}, titles))
                
                set(findall(gcf, '-property', 'fontsize'), 'fontsize', 18)
                set(findall(gcf, '-property', 'linewidth'), 'linewidth', 1.2)
                
                % save
                %                 folderName = sprintf('VSS2023/fig/RC/%s/corr_among_params/%s/all%d/', nameEnergySource, namesFeature{ifeature}, nLoc);
                %                 folderDir = dir(folderName); if isempty(folderDir), mkdir(folderName), end
                %                 saveas(gcf, sprintf('%sn%d_p%dp%d_%s_%s%d.jpg', folderName, nsubj, icol, irow, figTitle, text_sig, nsig))
                
                folderName = sprintf('%s/%s/corr_among_TunC/%s/%s/', nameFigFolder, name_numFilters_Fitting, nameFileLoc, namesFeature{ifeature});
                folderDir = dir(folderName); if isempty(folderDir), mkdir(folderName), end
                saveas(gcf, sprintf('%sn%d_p%dp%d_%s%d.jpg', folderName, nsubj, icol, irow, text_sig, nsig))
                
            end % if irow ~= icol
        end % icol
    end% irow
end % ii
%     close all

