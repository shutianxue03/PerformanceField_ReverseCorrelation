
perfThresh_plot = perfThresh_all(iPerf_plot)/100;

%% extract raw data
cst_log_unik = cst_log_unik_all{iLoc};
nCorr = nCorr_all{iLoc};
nData = nData_all{iLoc};

%% plot pC as a fxn of log cst
if ~exist('scaling', 'var'), scaling = 1; end
for icst_unik = 1: length(cst_log_unik)
    plot(cst_log_unik(icst_unik), nCorr(icst_unik)/nData(icst_unik), 'ok','MarkerSize',1+round(nData(icst_unik))/3*scaling,'Linewidth',1, 'HandleVisibility', 'off');
end
curveX_log = log10(fit.curveX).';

%% plot PMF and extract thresh
%CHANGE BACK AFTER PLOTTING PILOT DATA
legends_models = {};

LL_med = getCI(LL_all{iLoc}, 1, 1);
LL_med_delta = LL_med - max(LL_med);

for iModel = 1:nModels
    
    if iModel==4, x = 10.^cst_log_unik; X = 10.^curveX_log; % weibull fit to linear values
    else, x = cst_log_unik; X = curveX_log;
    end
    
    %% extract medians and CIs
    [yfit_med, yfit_lb, yfit_ub] = getCI(yfit_all{iLoc}(:, iModel, :), 1, 1);
    [thresh_med, thresh_lb, thresh_ub] = getCI(thresh_all{iLoc}(:, iModel), 1, 1);
    %%%%%%%%%%%%
    
    % plot the median of predicted pC
    plot(curveX_log, yfit_med,[ '-', colors_allM{iModel}]) % do NOT change curveX to X!!
    % plot the 68% CI of predicted pC
    patch([curveX_log; flip(curveX_log)], [yfit_lb(:); flip(yfit_ub(:))], colors_allM{iModel}, 'facealpha', .1, 'EdgeColor', 'none', 'HandleVisibility', 'off')
    % plot the median of thresh
    plot([-3, thresh_med], [perfThresh_plot, perfThresh_plot], '-', 'color', ones(1,3)/2, 'HandleVisibility', 'off')
    plot([-3, 0], [perfThresh_plot, perfThresh_plot], '-', 'color', ones(1,3)/2, 'HandleVisibility', 'off')
    plot([thresh_med, thresh_med], [0, perfThresh_plot], ['-', colors_allM{iModel}], 'HandleVisibility', 'off')
%     yline(.75, '-', 'color', ones(1,3)*.7, 'HandleVisibility', 'off'); % 75%, which should converge with threshold estimated from titration
    
    %CHANGE BACK AFTER PLOTTING PILOT DATA
    legends_models{iModel} = sprintf('%s (%.1f) %.1f%%', PMF_models{iModel}(1), LL_med_delta(iModel), 100*10^thresh_med);
end % im=1:nModels


%%

yticks([0, .5,  perfThresh_all/100, 1])
xticks_log = -3:.2:0;
xlim(xticks_log([1,end]))
xticks(xticks_log), xtickangle(90)
xticklabels(round(10.^xticks_log*100, 1))
ylim([0, 1])

yline(.5, 'color', ones(1,3)*.5);
set(gca, 'XGrid', 'on', 'YGrid', 'off') % grid on

%CHANGE BACK AFTER PLOTTING PILOT DATA
legend(legends_models, 'Location', 'south', 'NumColumns', 2)



