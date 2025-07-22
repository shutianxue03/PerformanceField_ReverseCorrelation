

%%
% compare the (log) RT at 5 loc between:
% 1. PRS vs. ABS
% 2. correct vs. wrong
% 3. say yes vs. no

%%
namesSAT = {'PRS', 'ABS'; 'Correct', 'Incorrect'; 'Resp=YES', 'Resp=NO'};
namesVar = {'trial type', 'Loc', 'Interaction'};
varInd_loc = [ones(1, nAllTrials); ones(1, nAllTrials)*2; ones(1, nAllTrials)*3; ones(1, nAllTrials)*4; ones(1, nAllTrials)*5];

figure('Position', [0 400 1200 400])
for n = 1:3
    subplot(1,3,n), hold on
    varInd_type = nan(nLoc5, nAllTrials);
    varList = varInd_type;
    
    for iLoc = 1:nLoc5
        if nsubj==1
            RT1_perLoc = RT_log_allT_perLoc{n, 1, iLoc};
            RT2_perLoc = RT_log_allT_perLoc{n, 2, iLoc};
        else
            RT1_perLoc = squeeze(RT_log_mean_allSubj(:, n, 1, iLoc));
            RT2_perLoc = squeeze(RT_log_mean_allSubj(:, n, 2, iLoc));
        end

        % ANOVA vector
        num1 = length(RT1_perLoc);
        num2 = length(RT2_perLoc);
        varInd_type(iLoc, :) = [ones(1, num1), ones(1, num2)*2];
        varList(iLoc, :) = [RT1_perLoc', RT2_perLoc'];
        
        % get mean and SEM
        RT1_mean = mean(RT1_perLoc);
        RT2_mean = mean(RT2_perLoc);
        RT1_SEM = std(RT1_perLoc)/sqrt(num1);
        RT2_SEM = std(RT2_perLoc)/sqrt(num2);
        
        RT_mean_perLoc(iLoc, [1,2]) = [RT1_mean, RT2_mean]; 
        RT_SEM_perLoc(iLoc, [1,2]) = [RT1_SEM, RT2_SEM]; 

        % plot
        errorbar(1:2, [RT1_mean, RT2_mean],[RT1_SEM, RT2_SEM], 'o-', 'color', colors_comb(iLoc, :), 'CapSize', 0)
    end
    
    ax =gca;
    ax.YGrid = 'on';
    ymin = min(RT_mean_perLoc(:)) - max(RT_SEM_perLoc(:))*2;
    ymax = max(RT_mean_perLoc(:)) + max(RT_SEM_perLoc(:))*2;
    yticks_ = linspace(ymin, ymax, 5);
    yticklabels_  = round(exp(yticks_)*1000);
    yticks(yticks_)
    yticklabels(yticklabels_)
    ylim([ymin, ymax])
    
    xlim([.5,2.5])
    xticks(1:2), xticklabels(namesSAT(n,:))
    ylabel('RT (ms)')
    
    % ANOVA
    [~, tbl] = anovan(varList(:), {varInd_type(:), varInd_loc(:)}, 'model','interaction','varnames',{'trial type', 'loc'}, 'display', 'off');
    titles = cell(1,3);
    for i = 1:3, titles{i} = sprintf('%s: F(%d,%d) = %.2f, p = %.3f', namesVar{i}, tbl{i+1,3}, tbl{5,3}, tbl{i+1,6}, tbl{i+1,7}); end
    title(titles)
end

set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)

make_sgtitle('Speed-accuracy tradeoff', subjName, nB, nsubj, nAllTrials)


