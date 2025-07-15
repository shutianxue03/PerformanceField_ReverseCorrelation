
clear ydiffErr

wd = 5;
wd_bar = 2;
sz_marker_ave = 30;
sz_marker_idvd = 12;
sz_ticks = 30;
sz_ylabel= 28;

%%%%%%%%%%%%%%%
% ORI
% if paramMode==1, ticks{1, 1} = 0:.12:.24; else, ticks{1, 1} = 0:.12:.24; end % gain/peak amp
ticks{1, 1} = 0:.1:.2; % peak amp
ticks{1, 2} =  0:45:90; % trough ori
ticks{1, 3} = -.04:.04:.04; % trough amp
ticks{1, 4} = 10:35:80; % band
ticks{1, 5} = [-.1, 0, .1]; % baseline

% SF
switch ifamily_perF(2)
    case 12
        % if paramMode==1, ticks{2, 2} = 0:.12:.24; else, ticks{2, 2} = 0:.06:.12; end % gain vs. peak amp
        % if flag_plotOctave, ticks{2, 3} = .4:1.2:2.8; else, ticks{2, 3} = 1:.5:3.5; end % octave vs. cpd
        ticks{2, 2} = 0:.08:.16; % peak amp 1
        ticks{2, 3} = .3:.4:1.1; % bandwidth 1
        ticks{2, 4} = [.678, 1.339, 2]; % peak SF2 (on log2 scale)
        ticks{2, 5} = 0:.06:.12; % peak amp 2
        ticks{2, 6} = .1:.5:1.1; % bandwidth 2
        %%%%%%%%%%%%%%%
    case 3
        ticks{2, 2} = 0:.08:.16; % peak amp
        ticks{2, 3} = .3:.4:1.1; % bandwidth
        ticks{2, 4} = [-.5, 0, .5]; % baseline
        ticks{2, 5} = 0:.06:.12; % trunc
end
xaxis = axis_tuning{ifeature};
nfilters = length(xaxis);
ifamily = ifamily_perF(ifeature);
namesParams = namesTunC_unit_perF{ifamily, paramMode};
nparams_full = length(namesParams);

% for estP, only plot gain
if paramMode == 1, if ifeature==1, pp_all = 1; else, pp_all = 2; end
else, pp_all = 1:nparams_full;
end

for iparam = pp_all
    
    figure('Position', [0 200 400 300]),hold on
    pp_allSubj = params_allSubj{ifeature}(:, :, :, itype, iparam);
    %     pp_allSubj = params_allSubj{ifeature}(:, :, :, itype, 1) - params_allSubj{ifeature}(:, :, :, itype, 3);
    
    % specifics for SF
    if ifeature == 2 & (find(iparam == [1,4])), pp_allSubj = log2(pp_allSubj);end
    
    %% get median/CI and ave/SEM
    [params_med, ~, ~, params_neg, params_pos] = getCI(pp_allSubj, 1, 2, 1, 1, CI_ratio); % CI for errorbar in the scatter plot
    [params_ave, ~, ~, params_sem] = getCI(params_med, 2, 1);
    fprintf('%s %s F>P: %d\n', namesFeature{ifeature},  namesParams{iparam}, sum(params_med(:,1)>params_med(:,2)))
    if (ifeature==2) && (iparam==3)
        params_ave = median(params_med);
    end
    
    %% plot bar
    for iLoc = 1:nLoc
        errorbar(iLoc, params_ave(iLoc), params_sem(iLoc), '.', ...
            'color', colors_comb_(iLoc, :), 'CapSize', 0, 'linewidth', wd)
        plot(iLoc, params_ave(iLoc), 'o', ...
            'MarkerEdgeColor', colors_comb_(iLoc, :), 'MarkerSize', sz_marker_ave, 'linewidth', wd)
    end % end of iiLoc
    
    %% idvd data
    if nLoc==2
        buffer = .2;
        for isubj = 1:nsubj
            plot([1+buffer, 2-buffer], params_med(isubj, :), [markers_allSubj{isubj}, '-'], ...
                'color',ones(1,3)*.7, 'markerfacecolor', 'w', 'markeredgecolor', ones(1,3)*.7, ...
                'markersize', sz_marker_idvd, 'linewidth', wd_bar)
        end
    else
        for isubj = 1:nsubj, plot(1:nLoc, params_med(isubj, :), [markers_allSubj{isubj}, '-'], 'color', ones(1,3)*.5, 'markersize', sz_marker_idv), end
    end
    
    %% compare with ref
    if nLoc==2
        x1 = params_med(:, 1);
        x2 = params_med(:, 2);
        
        % compare trough depth to 0
        ref=nan;
        if (ifeature==1) && (ifamily_perF(1) == 8) && (iparam == 3)
            ref = 0;
        end
        % compare SF peak to 2 cpd
        if (ifeature==2) && (ifamily_perF(2) == 12) && sum(iparam == [1,4])
            ref = 1;
        end
        if ~isnan(ref)
            yline(ref, 'color', ones(1,3)/2, 'linewidth', 2);
            fprintf('Comparing p%d to %d:\n', iparam, ref)
            [~, p, ~, stats] = ttest(x1-ref); fprintf('   %s: t=%.2f, p=%.3f (%d/%d)\n', namesLocComb_{1}, stats.tstat, p, sum(x1>2), nsubj);
            [~, p, ~, stats] = ttest(x2-ref); fprintf('   %s: t=%.2f, p=%.3f (%d/%d)\n', namesLocComb_{2}, stats.tstat, p, sum(x2>2), nsubj);
        end
    end
    
    %% draw sem of diff
    if nLoc==2
        diff_sem = std(x1-x2)/sqrt(nsubj);
        [~, p,~, stats] = ttest(x1, x2);
        cohenD = fxn_getES(x1, x2);
        flag_sig=''; if p<.05, flag_sig='sig'; end
        
        text_ANOVA = print_nANOVA({'Loc'}, [x1', x2'], {[ones(1, nsubj), ones(1, nsubj)*2]}, 0);
        
        ydiffErr = ticks{ifeature, iparam}(end)-(ticks{ifeature, iparam}(end) - ticks{ifeature, iparam}(1))/7;
        if (ifeature==1) && (iparam==2), ydiffErr = 10; end
        errorbar(1.5, ydiffErr, diff_sem, 'k', 'CapSize', 0, 'linewidth', wd)
        plot([1,2], [ydiffErr, ydiffErr], 'k-', 'linewidth', wd)
        string_s = getString_starts(p);
    end
    
    
    xticks([])
    buffer = .3;
    xlim([1-buffer, nLoc+buffer])
    
    %     if (~flag_plotOctave) && (ifeature == 2) && ~isempty((find(iparam == [1, 3]))),
    %         yticklabels(round(2.^ticks{ifeature, iparam}, 2));
    %     end
    if (ifeature==2) & (find(iparam==[1,4])), yticklabels(round(2.^ticks{ifeature, iparam}, 1)), end
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
    ax.LineWidth = wd;
    
    ylabel_ = sprintf('%s %s', namesFeature{ifeature}, namesParams{iparam});
    %     ylabel_ = sprintf('ORI peak-trough');
    
    if flag_plotStats
        title(sprintf('[p%d-%s]\nt(%d)=%.2f, p=%.3f, d=%.2f (%d/%d)', ...
            iparam, namesType{itype}, stats.df, stats.tstat, p, cohenD, sum(x1>=x2), nsubj), ...
            'fontsize', 10)
        %         ylabel(ylabel_, 'fontsize', sz_ylabel)
    end
    
    %     ylabel('distance to 2 cpd (log)') % delete
    % save
    folderName = sprintf('%s/%s/params/%s/%s/', nameFigFolder, name_numFilters_Fitting, nameFileLoc, namesType{itype});
    folderDir = dir(folderName); if isempty(folderDir), mkdir(folderName), end
    saveas(gcf, sprintf('%sn%d_%s%d_%s_%s.jpg', folderName, nsubj, namesFeature{ifeature}, iparam, namesParamsMode{paramMode}, flag_sig))
    %     saveas(gcf, sprintf('XueCarrasco/temp/params_L%d%d_n%d_%s%d.jpg', iLocComb_all,  nsubj, namesFeature{ifeature}, iparam)) % delete
    
end % end of iparam
