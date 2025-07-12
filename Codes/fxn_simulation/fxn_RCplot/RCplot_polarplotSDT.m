
figure('Position', [0 300 nLoc*200 200])

for ip = 1:ntitles % loop through each SDT measurements
    polarAxesHandle = subplot(1, ntitles, ip);
    ax = gca;
    ax.XTick = [];
    ax.YTick = [];
    polaraxes('Units',polarAxesHandle.Units,'Position',polarAxesHandle.Position)
    hold on
    if isubjPlot
        switch ip, case 1, x = dprime_allLoc; case 2, x = criterion_allLoc; case 3, x = RT_allLoc; end
        xmin = min(x)-nanstd(x);
        xmax = max(x)+nanstd(x);
        
        polarplot(polarAng, x(polarInd), 'color', colorSubj)
        polarplot(polarAngCircle, ones(1,1e3) * x(1), 'color', colorSubj)
        polarplot([polarAngFull;polarAngFull], repmat([xmin;xmax], 1,5), 'color', [.5,.5,.5])
        ax = gca;
        ax.ThetaTick = [];
        rlim([xmin, xmax])
    else
        switch ip, case 1, x = dprime_allSubj;case 2, x = criterion_allSubj; case 3, x = RT_allSubj; end
        % only ave data
        x_ave = nanmean(x,1); if ip == 3, x_ave = nanmedian(x,1);end
        xmin = min(x_ave)-nanstd(x_ave);
        xmax = max(x_ave)+nanstd(x_ave);
        polarplot(polarAngFull, x_ave(polarIndFull),'color', 'k')
        polarplot(polarAngCircle, ones(1,1e3) * x_ave(1), 'color', 'k')
        
        ax = gca;
        ax.ThetaTick = [];
        polarplot([polarAngFull;polarAngFull], repmat([xmin;xmax], 1,5), 'color', [.5,.5,.5])
        rlim([min(x_ave)-sum(abs(x_ave))/20, max(x_ave)+sum(abs(x_ave))/20])
    end
    title(SDT_titles{ip})
end

set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)