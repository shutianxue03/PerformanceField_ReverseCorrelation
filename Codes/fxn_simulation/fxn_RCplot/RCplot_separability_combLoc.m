

%% compare two locations
figure('Position', [0 200 ntypes*400 ncomb6*300])
for itype = 1:ntypes
    for icomb = 1:ncomb6
        subplot(ncomb6, ntypes, itype + ntypes*(icomb-1))
        hold on
        ax = gca; ax.YGrid = 'on';
        
        for iline = 1:2
            % get ave and CI/SEM
            if nsubj == 1
                if nB>1
                    [sep_ave, sep_lb, sep_ub, sep_SEM_neg, sep_SEM_pos] = getCI(seprblity(:, itype, combInd(icomb, iline)), errType, ntrialsProp);
                else
                    sep_ave = squeeze(seprblity(:, itype, combInd(icomb, iline)));
                    sep_SEM_neg = 0;
                    sep_SEM_pos = sep_SEM_neg;
                end
            else
                [sep_ave, ~, ~, sep_SEM_neg, sep_SEM_pos] = getCI(seprblity(:, itype, combInd(icomb, iline)), 2, 1, 1, ntrialsProp);
            end
            
            % plot group average
            bar(iline, sep_ave, 'FaceColor', 'w','EdgeColor', colors_comb(combInd(icomb, iline), :), 'BarWidth', .5, 'HandleVisibility', 'off')
            errorbar(iline, sep_ave, sep_SEM_neg, sep_SEM_pos, '.', 'color', colors_comb(combInd(icomb, iline), :), 'CapSize', 0, 'HandleVisibility', 'off')

        end % end of iline
        
        % plot IDVD data
        if nsubj>1
            for isubj_p = 1:nsubj
                plot(1:2, squeeze(seprblity(isubj_p, itype, [combInd(icomb, 1), combInd(icomb, 2)])), [marks_allSubj_full{isubj_p}, '-'], 'color', ones(1,3)*.5, 'linewidth', .5)
            end
        end
        
            % ANOVA
%             [~, tbl] = anova1(squeeze(seprblity(:, itype, :)), {}, 'off');
%             anova1_output = sprintf('F(%d, %d) = %.3f, p = %.3f', tbl{2,3}, tbl{3,3}, tbl{2,5}, tbl{2,6});
%             
            xlim([.2, 2.8])
            xticks(1:2)
            xticklabels({namesLocComb{combInd(icomb, :)}})
            
            ylim([0,1]), yticks(0:.2:1)
            if icomb == 1, title(namesType{itype}), end
            if icomb + itype == 2, ylabel('correlation coefficient'), end
%             if itype == 1 && nsubj>1, legend(subjList, 'Location', 'best','Orientation','horizontal' ,'NumColumns',2), end

    end % end of icomb
end % end of itypes
set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)

make_sgtitle([title_, ' - Separability'], subjName, nB, nsubj, nAllTrials)

