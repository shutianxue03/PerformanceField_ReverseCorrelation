
clear ydiffErr

wd_errbar1 = 4;
wd_errbar2 = 2;
wd_border = 2;
sz_marker1 = 20;
sz_marker2 = 12;
sz_title = 20;
%%%%%%%%%%%%%%%
RCplot_setting_margParams % define ticks,limits and the loc of texts
%%%%%%%%%%%%%%%

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
    
    figure('Position', [0 200 400 400]),hold on
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
    %     if
    [params_ave, ~, ~, params_sem] = getCI(params_med, 2, 1);
    fprintf('%s %s F>P: %d\n', namesFeature{ifeature},  namesParams{iparam}, sum(params_med(:,1)>params_med(:,2)))
    if (ifeature==2) && (iparam==3), 
        params_ave = median(params_med); end
    
    %% plot bar
    %     hold on%, box on
    
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
    
    %% idvd data
    if nLoc==2
        for isubj = 1:nsubj
            plot([1.3, 1.7], params_med(isubj, :), [markers_allSubj{isubj}, '-'], ...
                'color',ones(1,3)*.5, 'markerfacecolor', 'w', 'markeredgecolor', ones(1,3)*.5, ...
                'markersize', sz_marker2, 'linewidth', 1.2)
        end
    else
        for isubj = 1:nsubj, plot(1:nLoc, params_med(isubj, :), [markers_allSubj{isubj}, '-'], 'color', ones(1,3)*.5, 'markersize', sz_marker2), end
    end
    
    %% compare and draw sem of diff
    if nLoc==2
        x1 = params_med(:, 1);
        x2 = params_med(:, 2);
        diff_sem = std(x1-x2)/sqrt(nsubj);
        [~, p,~, stats] = ttest(x1, x2);
        ydiffErr = ticks{ifeature, iparam}(end)-(ticks{ifeature, iparam}(end) - ticks{ifeature, iparam}(1))/4;
         flag_sig=''; if p<.05, flag_sig='sig'; end
        text_ANOVA = print_nANOVA({'Loc'}, [x1', x2'], {[ones(1, nsubj), ones(1, nsubj)*2]}, 0);
        errorbar(1.5, ydiffErr, diff_sem, 'k', 'CapSize', 0, 'linewidth', wd_errbar1)
        plot([1,2], [ydiffErr, ydiffErr], 'k-', 'linewidth', wd_errbar1)
        string_s = getString_starts(p);
        
    end
    
    if (ifeature==2) && (iparam == 1), yline(1, 'color', ones(1,3)*.5, 'linewidth', 2); end
    
    axis square
    xticks(1:nLoc)
    xticklabels(namesLocComb(iLocComb_all))
    xlim([.3, nLoc+.7])

    if (~flag_plotOctave) && (ifeature == 2) && ~isempty((find(iparam == [1, 3]))),
        yticklabels(round(2.^ticks{ifeature, iparam}, 2));
    end
    if (ifeature==2) && (iparam==1), yticklabels(round(2.^ticks{ifeature, iparam}, 2)), end
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

    ax = gca;
    ax.XAxis.FontSize = sz_ticks;
    ax.YAxis.FontSize = sz_ticks;
    ax.LineWidth = wd_border;
    
    ax = gca;
    ax.XAxis.FontSize = sz_ticks;
    ax.YAxis.FontSize = sz_ticks;
    ax.LineWidth = wd_border;
    ylabel_ = sprintf('%s %s', namesFeature{ifeature}, namesParams{iparam});
    ylabel(ylabel_)
    set(findall(gcf, '-property', 'fontsize'), 'fontsize',30)
    %     set(findall(gcf, '-property', 'LineWidth'), 'LineWidth',2)
    
    
    %               title(sprintf('[%s]\nt(%d)=%.2f, p=%.3f\n%s\n%s', namesType{itype}, stats.df, stats.tstat, p, text_ANOVA, ylabel_), 'fontsize', sz_title)
    
    title(sprintf('[%s]\nt(%d)=%.2f, p=%.3f\n%s', namesType{itype}, stats.df, stats.tstat, p, text_ANOVA), 'fontsize', sz_title)
    %% save
    folderName = sprintf('%s/%s/params_simple/%s/%s/', nameFigFolder, nameEnergySource, nameFileLoc, namesType{itype});
    folderDir = dir(folderName); if isempty(folderDir), mkdir(folderName), end
    saveas(gcf, sprintf('%sn%d_%s%d_%s_%s.jpg', folderName, nsubj, namesFeature{ifeature}, iparam, namesParamsMode{paramMode}, flag_sig))
    
end % end of iparam
