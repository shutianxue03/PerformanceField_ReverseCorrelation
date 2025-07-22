
figure('Position', [0 200 1000 1000])
for itype = 1:ntypes
    for iLoc = 1:nLoc
        subplot(5, ntypes,(iLoc -1)*ntypes + itype), hold on
        plot(filterSF_all_log, energy_byPhase_allLoc{itype, iLoc})
        plot(filterSF_all_log, mean(energy_byPhase_allLoc{itype, iLoc}), 'k-', 'linewidth', 2)
        xline(log2(2), 'color', [.5,.5,.5]);
        yline(0, 'color', [.5,.5,.5]);
        xticks(log2([1,2,4])), xticklabels([1,2,4])
        xlim(log2([1,4])), ylim([-1,1])
        if iLoc == 1, title(['target-',typeNames{itype},' trials'])
            if itype == 3, leg = legend({'0', 'pi/4', 'pi/2' ,'3pi/4' ,'pi', 'ave'}, 'location', 'bestoutside'); title(leg,'add. phase')
            elseif itype == 1, xlabel('SF channels (cpd)'), ylabel('kernel (a.u.)')
            end
        end
    end
end
set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
sgtitle(['phase check - ', subjName, ' (each row is a loc)'], 'fontsize', 25)

