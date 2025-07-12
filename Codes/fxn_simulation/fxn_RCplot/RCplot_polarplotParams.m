figure('Position', [nLoc * 200 300 nLoc * 200 200])
for ip = 1:nparams
    polarAxesHandle = subplot(1, nparams, ip);
    ax = gca;
    ax.XTick = [];
    ax.YTick = [];
    polaraxes('Units',polarAxesHandle.Units,'Position',polarAxesHandle.Position)
    hold on
    if isubjPlot
        for itype = 1:ntypes
            x = squeeze(params(itype, :, ip));
%             xmin = min(x)-nanstd(x);
%             xmax = max(x)+nanstd(x);
            polarplot(polarAng, x(polarInd), 'color', colors{itype})
        end
%         polarplot(polarAngCircle, ones(1,1e3) * x(1), 'color', colorSubj)
%         polarplot([polarAngFull;polarAngFull], repmat([xmin;xmax], 1,5), 'color', [.5,.5,.5])
        ax = gca;
        ax.ThetaTick = [];
%         rlim([xmin, xmax])
    else
        % only the ave data
        x = squeeze(params(:,itype, :,ip));
        x_ave = nanmean(x,1);
        xmin = min(x_ave)-nanstd(x_ave);
        xmax = max(x_ave)+nanstd(x_ave);
        polarplot(polarAngFull, x_ave(polarIndFull),'color', 'k')
        polarplot(polarAngCircle, ones(1,1e3) * x_ave(1), 'color', 'k')
        polarplot([polarAngFull;polarAngFull], repmat([xmin;xmax], 1,5), 'color', [.5,.5,.5])
        ax = gca;
        ax.ThetaTick = [];
        rlim([xmin, xmax])
    end
    title(titlesParams{ip})
end

set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)