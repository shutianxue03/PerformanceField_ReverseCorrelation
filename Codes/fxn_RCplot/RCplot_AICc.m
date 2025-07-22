
function RCplot_AICc(nsubj, nmodels, paramInd_all, delta_allSubj2,plotIDVD, flagLongXTicks, title_)
sz_font = 8;

%% get mean and sem
delta_ave = mean(delta_allSubj2, 1);
delta_sem = std(delta_allSubj2, [], 1)/sqrt(nsubj);

%% get ranking
[delta_ave_s, sort_ind] = sort(delta_ave);
delta_sem_s = delta_sem(sort_ind);
% delta_allSubj2_s = delta_allSubj2(:, sort_ind);
paramInd_all_s = paramInd_all(sort_ind, :);

%% select good/best model by comparing AIC with 2
% if plotIDVD
%     [h,p] = ttest(delta_allSubj2_s, 2, 'tail', 'left');
% else
%     h=(delta_allSubj2_s <= 2);
% end

%
[a,b] = min(delta_ave_s);
hold on

% plot ave
for imodel = 1:nmodels
    if h(imodel) == 1, color_optimal = 'r';
    else
        if delta_ave_s(imodel) <= 2, color_optimal = 'm';
        elseif imodel == b, color_optimal = [1 .5 0];
        else, color_optimal = 'k';
        end
    end
    % bar+text point in red if smaller than 2
    bar(imodel, delta_ave_s(imodel), 'FaceColor', 'w', 'EdgeColor', color_optimal)
    text(imodel, delta_ave_s(imodel)+std(delta_allSubj2(:))/10, sprintf('%.1f', delta_ave_s(imodel)), 'HorizontalAlignment', 'center', 'fontsize' ,sz_font, 'color', color_optimal)
end

%% plot IDVD data
% if plotIDVD
%     errorbar(1:nmodels, delta_ave_s, delta_sem_s, 'k.', 'CapSize', 0)
%     for isubj = 1:nsubj
%         for imodel = 1:nmodels
%             % data point in red if smaller than 2
%             if delta_allSubj2_s(isubj, imodel) <= 2, color_optimal = 'r';
%             else, color_optimal = ones(1,3)/2;
%             end
%             plot(imodel+randn/8, delta_allSubj2_s(isubj, imodel), marks_allSubj_full{isubj}, 'color', color_optimal)
%         end
%     end
% end

%%
ylim([0,max(delta_allSubj2_s(:))])
% ylim([0,20])
yline(2, 'k--');

%% x tick labels
for imodel = 1:nmodels
    %     xticklabels_short{imodel} = sprintf('M%d', imodel);
    %     xticklabels_long_ = xticklabels_short{imodel};
    xticklabels_long_ = [];
    for ip = 1:nparams_full
        xticklabels_long_ = sprintf('%s\\newline%d', xticklabels_long_ , paramInd_all_s(imodel, ip));
    end
    xticklabels_long{imodel} = xticklabels_long_;
end

xticks(1:nmodels)
if flagLongXTicks, xticklabels(xticklabels_long)
else, xticklabels(xticklabels_short)
end

ax = gca;
ax.YGrid = 'ON'; 
ax.LineWidth=2;
ax.TickLength=[0 0];%



