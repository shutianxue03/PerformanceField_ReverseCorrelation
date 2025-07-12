
%% plot IDVD data of 5 locations
if nsubj>1
    
    figure('Position', [0 200 ntypes*400 300])
    for itype = 1:ntypes
        
        %     [sep_ave, ~, ~, sep_SEM_neg, sep_SEM_pos] = getCI(seprblity(:, itype, :), 2, 1, 1, ntrialsProp);
        
        subplot(1,ntypes, itype), hold on
        for iLoc = 1:nLoc5
            plot(1:nsubj, seprblity(:, itype, iLoc), '-o', 'color', colors_comb(iLoc, :))
            if nB>1
                %                 errorbar
            end
        end
        xlim([.5, nsubj+.5])
        xticks(1:nsubj)
        xticklabels(1:nsubj)
        ylabel('correlation coefficient')
        ylim([0,1]), yticks(0:.2:1)
        title(namesType{itype})
        
    end % end of itype
    set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
    make_sgtitle([title_, ' - Separability'], subjName, nB, nsubj, nAllTrials)
    
end


%% compare 5 locations
figure('Position', [0 200 ntypes*400 300])
for itype = 1:ntypes
    subplot(1, ntypes, itype)
    hold on
    ax = gca; ax.YGrid = 'on';
    
    % get ave and CI/SEM
    if nsubj == 1
        if nB > 1
            [sep_ave, sep_lb, sep_ub, sep_SEM_neg, sep_SEM_pos] = getCI(seprblity(:, itype, :), errType, ntrialsProp);
        else
            sep_ave = squeeze(seprblity(:, itype, :));
            sep_SEM_neg = zeros(1,nLoc8);
            sep_SEM_pos = zeros(1,nLoc8);
        end
    else
        [sep_ave, ~, ~, sep_SEM_neg, sep_SEM_pos] = getCI(seprblity(:, itype, :), 2, 1, 1, ntrialsProp);
    end
    
    % plot group average
    for iLoc = 1:nLoc5
        bar(iLoc, sep_ave(iLoc), 'FaceColor', 'w','EdgeColor', colors_comb(iLoc, :), 'BarWidth', .5, 'HandleVisibility', 'off')
        errorbar(iLoc, sep_ave(iLoc), sep_SEM_neg(iLoc), sep_SEM_pos(iLoc), '.', 'color', colors_comb(iLoc, :), 'CapSize', 0, 'HandleVisibility', 'off')
    end
    
    % plot IDVD data
    if nsubj>1
        for isubj_p = 1:nsubj
            plot(1:nLoc5, squeeze(seprblity(isubj_p, itype, 1:nLoc5)), [marks_allSubj_full{isubj_p}, '-'], 'color', ones(1,3)*.7, 'linewidth', .5)
        end
    end
    
    % ANOVA
    [~, tbl] = anova1(squeeze(seprblity(:, itype, :)), {}, 'off');
    anova1_output = sprintf('F(%d, %d) = %.3f, p = %.3f', tbl{2,3}, tbl{3,3}, tbl{2,5}, tbl{2,6});
    
    xlim([.5, nLoc5+.5])
    xticklabels(namesLoc2D)
    ylabel('correlation coefficient')
    ylim([0,1]), yticks(0:.2:1)
    title(namesType{itype})
    if itype == 1 && nsubj>1, legend(subjList, 'Location', 'best','Orientation','horizontal' ,'NumColumns',2), end
end

set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
make_sgtitle([title_, ' - Separability'], subjName, nB, nsubj, nAllTrials)

