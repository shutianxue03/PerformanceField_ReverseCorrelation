
%  check the corr among estP/tunC
% also look at corr between ORI and SF tunC

% SX_RC1_setting
clc
close all

labels_estP_ORI = namesTunC_unit_perF{ifamily_perF(1), 1};
labels_estP_SF = namesTunC_unit_perF{ifamily_perF(2), 1};
labels_tunC_ORI = namesTunC_unit_perF{ifamily_perF(1), 2};
labels_tunC_SF = namesTunC_unit_perF{ifamily_perF(2), 2};

sz_marker=20;


clear margParams
tunC = cat(ndims(tunC_ORI), tunC_ORI, tunC_SF); nTunC_ORI = length(labels_tunC_ORI);
labels = [labels_tunC_ORI, labels_tunC_SF];

nparams = length(labels);
nLoc = size(tunC, 3);
indLoc_corr = repmat(1:nLoc, nsubj, 1);
indSubj_corr = repmat((1:nsubj)', 1, nLoc);

iplots = 1:nparams*nparams; iplots = reshape(iplots, [nparams, nparams])';

for irow = 1:nparams
    if irow<=nTunC_ORI, nameFeature_row='ORI'; irow_title= irow; else, nameFeature_row='SF'; irow_title = irow-nTunC_ORI; end
    % param_y
    [y, ~, ~, y_neg, y_pos] = getCI(tunC(:, :, :, itype, irow), 1, 2);
    y_ = y-mean(y);
    
    for icol = 1:nparams
        if icol<=nTunC_ORI, nameFeature_col='ORI'; icol_title= icol; else, nameFeature_col='SF'; icol_title = icol-nTunC_ORI; end
        % param_x
        [x, ~, ~, x_neg, x_pos] = getCI(tunC(:, :, :, itype, icol), 1, 2);
        x_ = x-mean(x);
        
        % partial corr
        [r_partial, p_partial] = partialcorr([x(:), y(:)], indLoc_corr(:));
        r_partial = r_partial(2,1); p_partial = p_partial(2,1);
        titles = sprintf('\n%.2f (%.3f)', r_partial, p_partial);
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
                titles = [titles, sprintf('\n%s: %.2f (%.3f)', namesLocComb_{iLoc}, r_, p_)];
                if p_<.05, nsig = nsig+1; end
            end % iLoc
            
            % linear regression
            lm = polyfit(x_(:), y_(:), 1);
            x_lm2 = linspace(min(x_(:)) - std(x_(:)), max(x_(:)) + std(x_(:)), 2);
            yfit = polyval(lm, x_lm2);
            plot(x_lm2, yfit, 'k', 'linewidth', 2, 'handlevisibility', 'off');
            
            %
            axis square
            xlabel(sprintf('%s %s', nameFeature_col, labels{icol}))
            ylabel(sprintf('%s %s', nameFeature_row, labels{irow}))
            title(titles)
            
            set(findall(gcf, '-property', 'fontsize'), 'fontsize', 18)
            set(findall(gcf, '-property', 'linewidth'), 'linewidth', 1.2)
            
            % save
            if strcmp(nameFeature_col, 'ORI') && strcmp(nameFeature_row, 'ORI')
                CrossOrWithinFeature = 'ORI';
            elseif strcmp(nameFeature_col, 'SF') && strcmp(nameFeature_row, 'SF')
                CrossOrWithinFeature = 'SF';
            else
                CrossOrWithinFeature = 'ORI_vs_SF';
            end
            folderName = sprintf('%s/%s/corr_among_TunC/%s/%s/', nameFigFolder, name_numFilters_Fitting, nameFileLoc, CrossOrWithinFeature);
            folderDir = dir(folderName); if isempty(folderDir), mkdir(folderName), end
            saveas(gcf, sprintf('%sn%d_%s%d%s%d_%s%d.jpg', folderName, nsubj, nameFeature_col, icol_title, nameFeature_row, irow_title, text_sig, nsig))
            
        end % if irow ~= icol
    end % icol
end% irow
% end % ii
%     close all

