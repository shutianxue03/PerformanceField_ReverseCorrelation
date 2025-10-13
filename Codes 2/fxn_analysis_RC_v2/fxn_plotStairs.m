
%%
thresh_stair_single = nan(1, nLocSingle);

for iLoc = 1:nLoc
    endpoints = nan(1,nStairs);
    for istair=1:nStairs
        staircase_log = log10(ccc_allTrials(ccc_allTrials(:, 1)==iLoc & ccc_allTrials(:, 2)==istair, 3));
        if istair<=4, endpoints(istair)=staircase_log(end); end
    end
    thresh_stair_single(iLoc) = nanmean(endpoints);
end % iLoc

%% collapse HM
thresh_stair = [thresh_stair_single, mean(thresh_stair_single([2,4]))];

%% plot
figure('Position', [0 0 2e3 2e3]),

for iLoc = 1:nLoc
    subplot(3,3, iplots5(iLoc)), hold on, grid on
    
    for istair=1:nStairs
        if istair==1, color = 'r'; % staircase for left-tilted
        elseif istair==2, color = 'b'; % staircase for left-tilted
        end
        staircase_log = log10(ccc_allTrials(ccc_allTrials(:, 1)==iLoc & ccc_allTrials(:, 2)==istair, 3));
        
        plot(staircase_log, '.-', 'Color', color)
        if istair<=4, endpoints(istair)=staircase_log(end); end
        
    end % istair
    yline(mean(endpoints), 'k-', 'LineWidth',3);
    xlabel('trial#')
    ylabel('contrast (%)')
    yticks(-1:.25:0)
    yticklabels(round(10.^(-1:.25:0)*100, 1))
    ylim([-1, 0]);
    endpoints = 100*10.^(endpoints);
    title(sprintf('[L%d]  %.1f%%', iLoc, mean(endpoints)))
    
end % iLoc

set(findall(gcf, '-property', 'fontsize'), 'fontsize',12)
set(findall(gcf, '-property', 'LineWidth'), 'LineWidth',1.5)
sgtitle(subjName)
saveas(gcf, sprintf('%s/%s_stair.jpg', nameFolder_fig_thresh, subjName))
