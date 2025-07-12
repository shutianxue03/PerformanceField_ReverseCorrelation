
% scatter plot
ticks_allM = {-.1:.4:1.5,-.02:.06:.22, .5:.1:1};

%% figure setting
wd_border = 5;
sz_title = 20;
sz_label = 50; %40
sz_ticks = 50; % 30
sz_marker = 30;

for ip = 1:(nNoise-1)
    
    figure('Position', [0 0 1.2e3 1e3]), hold on, box on
    
    if flag_plotCI
        errorbar(NOM_params_med_allSubj(:, 1, ip), NOM_params_med_allSubj(:, 2, ip), ...
            NOM_params_neg_allSubj(:, 2, ip), NOM_params_pos_allSubj(:, 2, ip), ...
            NOM_params_neg_allSubj(:, 1, ip), NOM_params_pos_allSubj(:, 1, ip), ...
            '.', 'color', ones(1,3)*.5,'Capsize', 0, 'HandleVisibility','off')
    end
    for isubj = 1:nsubj
        plot(NOM_params_med_allSubj(isubj, 1, ip), NOM_params_med_allSubj(isubj, 2, ip), ...)
            markers_allSubj{isubj}, 'MarkerSize', sz_marker, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', ones(1,3)*.5, 'linewidth', wd_border)
    end
    nSubj_fp = sum(NOM_params_med_allSubj(:, 1, ip)<NOM_params_med_allSubj(:, 2, ip));
    % group average
    [ave_L1, ~, ~, sem_L1] = getCI(NOM_params_med_allSubj(:, 1, ip), 2, 1);
    [ave_L2, ~, ~, sem_L2] = getCI(NOM_params_med_allSubj(:, 2, ip), 2, 1);
    errorbar(ave_L1, ave_L2, sem_L1, 'horizontal', 'm.', 'Capsize', 0, 'linewidth', wd_border*2)
    errorbar(ave_L1, ave_L2, sem_L2, 'vertical', 'm.', 'Capsize', 0, 'linewidth', wd_border*2)
    
    % print range and number of observers whose fovea < peri
    fprintf('Number of subj whose noise is lower at Fovea: %d \n%s\n    %s %.4f to %.4f (%.2f +- %.2f)\n    %s %.4f to %.4f (%.2f +- %.2f)\n\n', ...
        nSubj_fp, namesNOMparams{ip}, ...
        namesLocComb{iLocComb_all(1)}, min(NOM_params_med_allSubj(:, 1, ip)), max(NOM_params_med_allSubj(:, 1, ip)), ave_L1, sem_L1, ...
        namesLocComb{iLocComb_all(2)}, min(NOM_params_med_allSubj(:, 2, ip)), max(NOM_params_med_allSubj(:, 2, ip)), ave_L2, sem_L2)

    %         if ip == 1, legend(subjList, 'Location', 'southeast', 'NumColumns', 3), end
    plot(ticks_allM{ip}([1, end]), ticks_allM{ip}([1, end]), 'k-', 'handlevisibility', 'off', 'linewidth', wd_border)
    xlim(ticks_allM{ip}([1, end]))
    ylim(ticks_allM{ip}([1, end]))
    xticks(ticks_allM{ip})
    yticks(ticks_allM{ip})
    xlabel(namesLocComb{iLocComb_all(1)})
    ylabel(namesLocComb{iLocComb_all(2)})
    axis square
    
    [~, p, CI, stats] = ttest(NOM_params_med_allSubj(:, 1, ip), NOM_params_med_allSubj(:, 2, ip));
    
    ax = gca;
    ax.XAxis.FontSize = sz_ticks;
    ax.YAxis.FontSize = sz_ticks;
    ax.LineWidth = wd_border;
    title(sprintf('%s\nt(%d) = %.3f, p = %.3f', namesParamsNOM{ip}, stats.df, stats.tstat, p), 'FontSize', sz_title)
    
    %% save
    folderName = sprintf('%s/temp%d/comp_params/%s/', nameFigFolder, templateType, nameFileLoc);
    folderDir = dir(folderName); if isempty(folderDir), mkdir(folderName), end
    saveas(gcf, sprintf('%sn%d_P%d_A%dB%d.jpg', folderName, nsubj, ip, iModelA, iModelB))
end % ip
