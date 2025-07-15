


ylabels = {'Hit rate', 'False alarm rate'};
titles = {'Signal-present trials', 'Signal-absent trials', 'All trials'};
xtext = 2;
wd_border=1; sz_ticks = 20; sz_title=20; sz_marker = 30/nLoc;

%% 2 way ANOVA (Loc x ebin on p(YES) for PRS/ABS separately)
ind_Loc2 = repmat(1:nLoc, [nsubj, 1, nbins_e]); % 2 stands for 2-way ANOVA
ind_ebin2 = [];
for ib = 1:nbins_e
    ind_ebin_ = ones(nsubj, nLoc)*ib;
    ind_ebin2 = cat(3, ind_ebin2, ind_ebin_);
end

%% 3-way ANOVA (iPrs x Loc x ebin on p(YES))
r_anova3 = nan(nsubj, nLoc, nbins_e, 2); %12x2x5x2
ind_iPrs3 = [];
for ip = 1:2
    ind_iPrs_ = ones(nsubj, nLoc, nbins_e) * ip;
    ind_iPrs3 = cat(4, ind_iPrs3, ind_iPrs_);
end
ind_Loc3 = nan([size(ind_Loc2), 2]);
ind_ebin3 = nan([size(ind_ebin2), 2]);

%%
[pYES_med, ~, ~, pYES_neg, pYES_pos] = getCI(pYES_tgt_allSubj_, 1, 2); % nsubj x ntypes x nLoc x nbins_e
[ebin_med, ~, ~, ebin_neg, ebin_pos] = getCI(ebin_tgt_allSubj_, 1, 2);

figure('Position', [0 300 1000 500]), hold on, box on

for itype_ = 1:2
    subplot(1,2,itype_), hold on, box on
    if itype_ == 1
        if flag_PatchMode==1, xticks_ = .08:.08:.32; xlim_ = [.064, .336]; % energy derived from target oatches
        else, xticks_ = 0:.04:.12; xlim_ = [-.008, .128]; % energy derived from noise patches
        end
        ytext = .12;
    else
        xticks_ = 0:.04:.12; xlim_ = [-.008, .128];
        ytext = .8;
    end
    
    [ebin_ave, ~, ~, ebin_sem] = getCI(ebin_med(:, :, itype_, :), 2,1);
    [pYES_ave, ~, ~, pYES_sem] = getCI(pYES_med(:, :, itype_, :), 2,1);
    
    for iLoc = 1:nLoc
        x = ebin_ave(iLoc, :);
        y = pYES_ave(iLoc, :);
        color = colors_comb_(iLoc, :);
        % raw data
        errorbar(x, y, ebin_sem(iLoc, :), 'horizontal','.', 'color', color, 'linewidth', wd_border, 'CapSize', 0, 'handlevisibility', 'off')
        errorbar(x, y, pYES_sem(iLoc, :), 'vertical',  '.', 'color', color, 'linewidth', wd_border, 'CapSize', 0, 'handlevisibility', 'off')
        plot(x, y, 'o', 'MarkerFaceColor', color, 'MarkerEdgeColor', 'w', 'MarkerSize', sz_marker, 'linewidth', wd_border)
        
        [pearsonr,p] = corr(x',y');
        ss = getString_starts(p);
        x_lm = linspace(xticks_(1), xticks_(end), 2);
        [beta, stats] = polyfit(x, y, 1);
        slope = beta(1);
        yfit = polyval(beta, x_lm);
        plot(x_lm, yfit, '.-', 'color', color, 'handlevisibility', 'off', 'linewidth', wd_border)
        
    end % iiLoc
    
    % ANOVA2 (3 bins x 5 loc on rFA and rHit separately)
    r_anova2 = squeeze(pYES_med(:, :, itype_, :)); % nsubj x nLoc x nbins_e
    
    % save as csv and run in JASP
    %     r_csv = [r_(:), ind_Loc2(:), ind_ebin2(:)];
    %     csvwrite('Data/r.csv', [1:3;r_csv])
    
    ylim([0 1]), yticks(0:.25:1)
    xlabel('Binned energy', 'FontSize', sz_ticks)
    ylabel(ylabels{itype_}, 'FontSize', sz_ticks)
    xticks(xticks_)
    xlim(xlim_)
    
    ax = gca;
    ax.XAxis.FontSize = sz_ticks;
    ax.YAxis.FontSize = sz_ticks;
    ax.LineWidth = wd_border;
    
    anova_text = print_nANOVA({'Loc', 'Energy'}, r_anova2(:), {ind_Loc2(:), ind_ebin2(:)}, nsubj);
    fprintf('%s\n%s\n', namesType{itype_}, anova_text)
    r_anova3(:, :, :,itype_) = r_anova2;
    ind_Loc3(:, :, :,itype_) = ind_Loc2;
    ind_ebin3(:, :, :,itype_) = ind_ebin2;
    
    title(anova_text)
end % end of itype_

% 3-way ANOVA (defined above)
varNames = {'iPRS', 'Loc', 'ebin'};
anova_text = print_nANOVA(varNames, r_anova3(:),  {ind_iPrs3(:), ind_Loc3(:), ind_ebin3(:)}, nsubj);

sgtitle(anova_text)

folderName = sprintf('%s/%s/ebin/', nameFigFolder, nameEnergySource);
folderDir = dir(folderName); if isempty(folderDir), mkdir(folderName), end

saveas(gcf, sprintf('%sn%d_%s.jpg', folderName, nsubj, nameFileLoc))


