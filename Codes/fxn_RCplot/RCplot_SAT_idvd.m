


%% conduct ANOVA
% nVars = ndims(RT);
[ntrials, nLoc] = size(RT_log); % ntrials (of all sessions) x nLoc

namesVar = {'trial type', 'Loc', 'Interaction'};

RT_ = RT_log(:);
corr_ = correct(:);
var2ind = repmat(1:nLoc,ntrials,1); var2ind = var2ind(:);
[~, tbl] = anovan(exp(RT_), {corr_, var2ind},'model','interaction','varnames',{'trial type', 'loc'}, 'display', 'off');

titles = cell(1,3);
for n = 1:3, titles{n} = sprintf('%s: F(%d,%d) = %.2f, p = %.3f', namesVar{n}, tbl{n+1,3}, tbl{5,3}, tbl{n+1,6}, tbl{n+1,7}); end

%% get mean and SD
RT_corr_mean = nan(nLoc, 1);
RT_inc_mean = RT_corr_mean;
RT_corr_sem = RT_corr_mean;
RT_inc_sem = RT_corr_mean;

for iLoc = 1:nLoc
    RT_corr = RT_log(boolean(correct(:,iLoc)), iLoc);
    RT_inc = RT_log(boolean(1-correct(:,iLoc)), iLoc);
    RT_corr_mean(iLoc) = mean(RT_corr);
    RT_inc_mean(iLoc) = mean(RT_inc);
    RT_corr_sem(iLoc) = std(RT_corr)/sqrt(ntrials);
    RT_inc_sem(iLoc) = std(RT_inc)/sqrt(ntrials);
end

%% plot
ymin = min([RT_corr_mean; RT_inc_mean]);
ymax = max([RT_corr_mean; RT_inc_mean]);

figure('Position', [200 0 400 400]), hold on
for iLoc = 1:nLoc
    errorbar(1:2, [RT_corr_mean(iLoc), RT_inc_mean(iLoc)], [RT_corr_sem(iLoc), RT_inc_sem(iLoc)],'o-', 'Capsize', 0, 'color', colors5Loc(iLoc, :))
end

yticks(linspace(ymin, ymax, 5))
yticklabels(round(exp(linspace(ymin, ymax, 5))*1000))
xticks([1,2]), xticklabels({'correct trials', 'incorrect trials'}), xtickangle(45)
ylabel('RT (ms)')
xlim([0,2.5])
legend(namesLoc2D, 'location', 'northwest')
title(titles)




