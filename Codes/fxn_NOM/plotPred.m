function plotPred(iplot, nsubj, iLoc_all, data, pred, title_, ylim_)
% data structure
%       data: nLoc x 1 (IDVD) // nsubj x nLoc (Group)
%       pred: ni//nsubj x nLoc
%       params_est: ni//nsubj x nLoc x nparams

nLoc = length(iLoc_all);
namesLocComb = {'Fovea', 'LHM', 'UVM', 'RHM', 'LVM', 'HM', 'VM', 'Peri'};

if size(data, 1) == nsubj; modeGroup = 1; else, modeGroup = 0; end

subplot(4,3, iplot), hold on

buffer = .2;
if modeGroup
    data_mean = mean(data, 1);
    data_sem = std(data, [], 1)/sqrt(nsubj);
    pred_mean = mean(pred, 1);
    pred_sem = std(pred, [], 1)/sqrt(nsubj);
    
    bar(1:nLoc, data_mean, 'FaceColor', 'none', 'EdgeColor', 'k')
    for isubj = 1:nsubj
        plot(1:nLoc, data(isubj, :), '-', 'color', ones(1,3)*.5)
        plot((1:nLoc) - buffer, pred(isubj, :), '-', 'color', [1,0,0]*.5)
    end
    errorbar(1:nLoc, data_mean, data_sem, 'k.', 'CapSize', 0)
    errorbar((1:nLoc)-buffer, pred_mean, pred_sem, 'ko', 'CapSize', 0)
else
    
    bar(1:nLoc, data)
    for iLoc = 1:nLoc
        pred_med = median(pred(:, iLoc), 1); % 1st col: Loc1; 2nd col: loc 2
        pred_CI68 = quantile(pred(:, iLoc), [.16, .84], 1);
        errorbar(iLoc, pred_med, pred_med - pred_CI68(1), pred_CI68(2) - pred_med, 'ko', 'CapSize', 0)
    end
end

xticks(1:nLoc), xticklabels(namesLocComb(iLoc_all))
xlim([0, nLoc+1])
% if nargin == 7, ylim(ylim_), end
title(title_)
end