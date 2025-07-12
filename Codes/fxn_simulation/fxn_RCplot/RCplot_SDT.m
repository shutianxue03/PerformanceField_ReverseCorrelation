
ntitles = length(namesSDT);
yticks_ = [.5, 1, 1.5; -1, 0, 1; 0, .15, .3];

figure('Position', [0 200 nLoc*200 300])
for ititle = 1:ntitles % loop through each SDT measurements
    subplot(1, ntitles, ititle), hold on, box on
    
    if plotAve
        switch ititle, case 1, x = dprime_allSubj; case 2, x = criterion_allSubj; case 3, x = exp(RT_log_allSubj);
            % plot IDVD data
            for isubj_p = 1:nsubj, plot((1:nLoc) +.1, x(isubj_p, :), [marks_allSubj{isubj_p}, '-'], 'color', ones(1,3)*.5), end
            % get ave and SEM
            x_ave = squeeze(nanmean(x,1));
            x_sem = squeeze(nanstd(x, [], 1))/sqrt(nsubj);
        end
    else
        switch ititle, case 1, x_ave = dprime_perLoc; case 2, x_ave = criterion_perLoc; case 3, x_ave = exp(RT_log_perLoc); end
        x_sem = zeros(1,nLoc);
    end
    
    % plot average
    for iLoc = 1:nLoc
        bar(iLoc, x_ave(iLoc), 'FaceColor', colors_comb(iLoc, :), 'facealpha', .3, 'handlevisibility', 'off', 'EdgeColor', 'w')
        errorbar(iLoc, x_ave(iLoc), x_sem(iLoc), '.', 'color', colors_comb(iLoc, :), 'markersize', 10, 'linewidth', 2)
    end
    
    if plotAve
        % stats
        [~, tbl] = anova1(x, {}, 'off');
        anova1_output = sprintf('F(%d, %d) = %.3f, p = %.3f', tbl{2,3}, tbl{3,3}, tbl{2,5}, tbl{2,6});
    else, anova1_output = '';
    end
    % extra lines
    if ititle == 1, yline(dprime_theo, 'k');  end
    if ititle == 2, yline(0, 'k'); end
    
    ylim(yticks_(ititle, [1,end]))
    yticks(yticks_(ititle, :))
    xticks(1:nLoc), xticklabels(namesLoc2D)
    ax = gca; set(ax, 'FontSize', 20)
    xlim([.5, nLoc+.5])
    if ititle == 2, legend(subjList, 'Location', 'best','Orientation','horizontal'), end
    
    title(sprintf('%s\n%s', namesSDT{ititle}, anova1_output), 'FontSize', 15)
end
if nsubj==1, sgtitle([subjName, '- SDT']), end

set(findall(gcf, '-property', 'FontSize'), 'FontSize', 15)
% set(findall(gcf, '-property', 'linewidth'), 'linewidth',2)

if nsubj>1, saveas(gcf, sprintf('publishedPDFs/fig/n%d_SDT.jpg', nsubj)), end

