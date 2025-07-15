
clear ydiffErr

wd_errbar1 = 4;
wd_errbar2 = 2;
wd_border = 2;
RCplot_setting_margParams % define ticks,limits and the loc of texts
sz_marker1 = 18;
sz_marker2 = 15;
sz_title = 20;

xaxis = axis_tuning{ifeature};
nfilters = length(xaxis);
xaxis_itp = axis_itp_tuning{ifeature};
ifamily = ifamily_perF(ifeature);
namesParams = namesParams_unit_perF{paramMode, ifeature};
nparams_full = length(namesParams);

% for estP, only plot gain
if paramMode == 1, if ifeature==1, pp_all = 1; else, pp_all = 2; end
else, pp_all = 1:nparams_full;
end

for iparam = pp_all
    
    figure('Position', [0 200 900 600]);
    pp_allSubj = params_allSubj{ifeature}(:, :, :, itype, iparam);
    
    % specifics for SF
    if ifeature == 2
        % SF peak, plot log scale
        if iparam == 1, pp_allSubj = log2(pp_allSubj);end
        % bandwidth, plot log OR linear scale (saved as the last element)
        if (iparam == 3) && ~flag_plotOctave, pp_allSubj = params_allSubj{ifeature}(:, :, :, itype, end); end
    end
    
    %% get median/CI and ave/SEM
    [params_med, ~, ~, params_neg, params_pos] = getCI(pp_allSubj, 1, 2, 1, 1, CI_ratio); % CI for errorbar in the scatter plot
    [params_ave, ~, ~, params_sem] = getCI(params_med, 2, 1);
    
    %% plot bar
    subplot(1,2,1), hold on, box on
    
    % extra line
    if (ifeature==1) && (iparam == 3)
        yline(0, 'color', ones(1,3)*.5, 'linewidth', 2);
    elseif (ifeature==2) && (iparam == 1)
        yline(1, 'color', ones(1,3)*.5, 'linewidth', 2);
    end
    
    for iLoc = 1:nLoc
        errorbar(iLoc, params_ave(iLoc), params_sem(iLoc), '.', ...
            'color', colors_comb_(iLoc, :), 'CapSize', 0, 'linewidth', wd_errbar1)
        plot(iLoc, params_ave(iLoc), 'o', 'MarkerFaceColor' , colors_comb_(iLoc, :), ...
            'MarkerEdgeColor', 'w', 'MarkerSize', sz_marker1, 'linewidth', wd_errbar2)
    end % end of iiLoc
    
    % idvd data
    if nLoc==2
        for isubj = 1:nsubj, plot([1.3, 1.7], params_med(isubj, :), [markers_allSubj{isubj}, '-'], 'color', ones(1,3)*.5), end
    else
        for isubj = 1:nsubj, plot(1:nLoc, params_med(isubj, :), [markers_allSubj{isubj}, '-'], 'color', ones(1,3)*.5), end
    end
    % compare and draw sem of diff
    if nLoc==2
        x1 = params_med(:, 1);
        x2 = params_med(:, 2);
        diff_sem = std(x1-x2)/sqrt(nsubj);
        [~, p,~, stats] = ttest(x1, x2);
        ydiffErr = ticks{ifeature, iparam}(end)-(ticks{ifeature, iparam}(end) - ticks{ifeature, iparam}(1))/4;
        
        text_ANOVA = print_nANOVA({'Loc'}, [x1', x2'], {[ones(1, nsubj), ones(1, nsubj)*2]}, 0);
        errorbar(1.5, ydiffErr, diff_sem, 'k', 'CapSize', 0, 'linewidth', wd_errbar1)
        plot([1,2], [ydiffErr, ydiffErr], 'k-', 'linewidth', wd_errbar1)
        string_s = getString_starts(p);
        if flag_plotStats
            %             text(1.5, ydiffErr+diff_sem*2, ...
            %                 sprintf('t(%d)=%.2f, p=%.3f\n%s', stats.df, stats.tstat, p, text_ANOVA), 'HorizontalAlignment', 'center', 'fontsize', sz_title)
            title(sprintf('[%s]\nt(%d)=%.2f, p=%.3f\n%s', namesType{itype}, stats.df, stats.tstat, p, text_ANOVA), 'fontsize', sz_title)
        end
    end
    
    if (ifeature==2) && (iparam == 1), yline(1, 'color', ones(1,3)*.5, 'linewidth', 2); end
    
    axis square
    xticks(1:nLoc)
    xlim([.3, nLoc+.7])
    xticklabels(xTLs)
    
    ylabel(sprintf('%s %s', namesFeature{ifeature}, namesParams{iparam}))
    if (~flag_plotOctave) && (ifeature == 2) && ~isempty((find(iparam == [1, 3]))), yticklabels(round(2.^ticks{ifeature, iparam}, 2));end
    yticks(ticks{ifeature, iparam})
    ylim(ticks{ifeature, iparam}([1, end]))
    
    xlabel('') % to ensure the same shape of the bar plot on the left & the scatter plot on the right
    
    % compare some estP/TunC to reference
    flag_compRef=0;
    if ifeature==1
        if iparam == 3 % ORI baseline/bottom
            flag_compRef = 1; ref = 0; tail = 'left';
        end
    else
        if iparam == 1 % peak SF
            flag_compRef = 1; ref = 1; tail = 'both';
        elseif iparam == 4 % SF baseline/bottom
            flag_compRef = 1; ref = 0; tail = 'left';
        end
    end
    
    if flag_compRef
        fprintf('%s %s\n', namesFeature{ifeature}, namesParams{iparam})
        for iLoc = 1:nLoc
            [~, p, ~, stats] = ttest(params_med(:, iLoc), ref, 'tail', tail);
            fprintf('    %s: t(%d) = %.2f, p = %.2f\n', namesLocComb_{iLoc}, stats.df, stats.tstat, p)
        end
    end
    ax = gca;
    ax.XAxis.FontSize = sz_ticks;
    ax.YAxis.FontSize = sz_ticks;
    ax.LineWidth = wd_border;
    
    
    %% plot scatter
    if nLoc == 2
        subplot(1,2,2), hold on, box on
        
        % extra line
        if (ifeature==1) && (iparam == 3)
            xline(0, 'color', ones(1,3)*.5, 'linewidth', 2);
            yline(0, 'color', ones(1,3)*.5, 'linewidth', 2);
        elseif (ifeature==2) && (iparam == 1)
            xline(1, 'color', ones(1,3)*.5, 'linewidth', 2);
            yline(1, 'color', ones(1,3)*.5, 'linewidth', 2);
        end
        
        % diagonal line
        plot(ticks{ifeature, iparam}([1,end]), ticks{ifeature, iparam}([1,end]), 'k-', 'linewidth', wd_errbar2)
        
        % IDVD data
        errorbar(params_med(:,1), params_med(:,2), ...
            params_neg(:,2), params_pos(:,2), ...
            params_neg(:,1), params_pos(:,1), ...
            '.', 'color', ones(1,3)*.5, 'CapSize', 0, 'linewidth', wd_errbar2)
        for isubj = 1:nsubj
            plot(params_med(isubj,1), params_med(isubj,2), markers_allSubj{isubj}, ...
                'linestyle', 'none', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', ones(1,3)*.5, 'linewidth', wd_errbar2, 'MarkerSize', sz_marker2)
        end
        
        % group ave (purple cross)
        errorbar(params_ave(1), params_ave(2), params_sem(1), 'horizontal', 'color', [1,0,1], 'linewidth', wd_errbar1, 'CapSize', 0)
        errorbar(params_ave(1), params_ave(2), params_sem(2), 'vertical', 'color', [1,0,1], 'linewidth', wd_errbar1, 'CapSize', 0)
        
        xticks(ticks{ifeature, iparam})
        yticks(ticks{ifeature, iparam})
        
        if (ifeature == 2) && (iparam == 1)
            xticklabels(round(2.^ticks{ifeature, iparam}, 2));
            yticklabels(round(2.^ticks{ifeature, iparam}, 2));
        end
        
        ylim(ticks{ifeature, iparam}([1, end]))
        xlim(ticks{ifeature, iparam}([1, end]))
        
        xlabel(xTLs{1})
        ylabel(xTLs{2})
        axis square
        
    end % if nLoc == 2
    
    ax = gca;
    ax.XAxis.FontSize = sz_ticks;
    ax.YAxis.FontSize = sz_ticks;
    ax.LineWidth = wd_border;
    set(findall(gcf, '-property', 'fontsize'), 'fontsize',25)
    
    %% save
    folderName = sprintf('%s/%s/params/%s/%s/', nameFigFolder, nameEnergySource, nameFileLoc, namesType{itype});
    folderDir = dir(folderName); if isempty(folderDir), mkdir(folderName), end
    saveas(gcf, sprintf('%sn%d_%s%d_%s.jpg', folderName, nsubj, namesFeature{ifeature}, iparam, namesParamsMode{paramMode}))
    
end % end of iparam
