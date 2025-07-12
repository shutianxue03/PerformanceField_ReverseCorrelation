for ifilterOri = ceil(nfiltersOri/2)-1:ceil(nfiltersOri/2)+1
    for ifilterSF = ceil(nfiltersSF/2)-1:ceil(nfiltersSF/2)+1
        figure('Position', [0 200 ntypes*400 800])
        for itype = 1:ntypes
            switch itype
                case 1, itrialStart = 1; itrialEnd = ntrials;
                case 2, itrialStart = ntrials+1; itrialEnd = ntrials*2;
                case 3, itrialStart = 1; itrialEnd = ntrials*2;
            end
           
            % plot resp rate vs. binned energy
            subplot(2, ntypes, itype), hold on, grid on
            slope_legends = cell(nLoc, 1);
            for iLoc = 1:nLoc
                [pC_allBins, energy_allBins] = checkEnerBins(energy_perLoc(iLoc, itrialStart:itrialEnd, ifilterOri, ifilterSF), data_both_B(itrialStart:itrialEnd, iLoc, 2), nbins_e);
                [beta, stats] = polyfit(energy_allBins, pC_allBins, 1); % linear fitting
                slope = beta(1); % beta = [slope, intercept]
                yfit = polyval(beta, energy_allBins);
                R2 = getR2(pC_allBins, yfit);
                slope_legends{iLoc} = sprintf('%s %.2f (%.1f%%)', namesLoc2D{iLoc}, slope, R2*100);
                plot(energy_allBins*100, pC_allBins, 'o', 'MarkerFaceColor',  colors5Loc(iLoc, :),  'MarkerEdgeColor', 'w', 'handlevisibility', 'off')
                plot(energy_allBins*100, yfit, '-', 'color' ,  colors5Loc(iLoc, :))
                slopes(iLoc) = slope;
                
            end
            legend(slope_legends, 'Location', 'best')
            ylim([0, 1])
            
            % plot the slope
            subplot(2, ntypes, itype+ntypes), hold on, grid on
            
            for iLoc = 1:nLoc
                bar(iLoc, slopes(iLoc), 'FaceColor', colors_comb(iLoc, :), 'FaceAlpha', .5)
            end
            ylim([-1,3])
            ylabel('slope')
            xticks(1:5), xticklabels(namesLoc2D), xtickangle(45)
            
        end
        set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
        sgtitle(sprintf('ORI=%d SF=%.2f', round(filtersOri_all(ifilterOri)-90), filtersSF_all(ifilterSF)), 'fontsize', 20)
    end
end