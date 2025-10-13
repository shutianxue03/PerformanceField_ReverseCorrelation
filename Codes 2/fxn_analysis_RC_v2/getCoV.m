clc, close all
itype=2;
iModelA = 1;
iLocComb_all = [6,5,3];nLoc=3;

nData = 9;
t_all_short = {'CS', 'pA', 'ORI peakAmp', 'ORI bandwidth', 'Pref SF', 'SF peakAmp', 'SF bandwidth', 'IN-Induced', 'IN-Const'};
t_all = {'Contrast sensitivity', 'Response consistency', 'ORI peak amp.', 'ORI bandwidth', 'Prefered SF', 'SF peak amp.', 'SF bandwidth', 'Induced noise', 'Constant noise'};
x_ticks_all = [-.5, .5; -.5, .5; .3, .6; .2, 1; .1, .5; .25, .5; 0, .8; .3, .8;.4, 1.4];
CoV_med = nan(nData, nLoc);
CoV_SEM_neg = CoV_med;
CoV_SEM_pos = CoV_med;
CoV_allB = cell(nData, nLoc);

for iData = 1:nData
    switch iData
        case 1, data = cs_allSubj;
        case 2, data = pA_allSubj;
        case 3, data = margTuningC_ORI_3allSubj(:, :, :, itype, 2);
        case 4, data = margTuningC_ORI_3allSubj(:, :, :, itype, 5);
        case 5, data = margTuningC_SF_3allSubj(:, :, :, itype, 1);
        case 6, data = margTuningC_SF_3allSubj(:, :, :, itype, 2);
        case 7, data = margTuningC_SF_3allSubj(:, :, :, itype, 3);
        case 8, data = NOM_params_est_allSubj(:, :, iModelA, :, :, 1);
        case 9, data = NOM_params_est_allSubj(:, :, iModelA, :, :, 2);
    end
    t = t_all_short{iData};
    data=squeeze(data);
    if iData>7, [nsubj, nLoc, nB] = size(data);
    else, [nsubj, nB, nLoc] = size(data);
    end
    
    fprintf('\n%s:\n    ', t)
    for iLoc = 1:nLoc
        if iData>7, data_ = squeeze(data(:, iLoc, :));
        else, data_ = data(:, :, iLoc);
        end
        CoV = std(data_, [], 1)./mean(data_, 1);
        
        [med, ~, ~, sem_neg, sem_pos] = getCI(CoV, 1, 2);
        CoV_med(iData, iLoc) = med;
        CoV_SEM_neg(iData, iLoc) = sem_neg;
        CoV_SEM_pos(iData, iLoc) = sem_pos;
        CoV_allB{iData, iLoc} = CoV;
        fprintf('L%d: %.2f [%.2f, %.2f];    ', iLocComb_all(iLoc), med, sem_neg, sem_pos)
    end % iLoc
end %iData

%%
clc, close all

buffer = .2;
nBins = 50;

for iData = 3:nData
    figure('Position', [0 0 400 300]), hold on
    
    for iLoc = 1:nLoc
        color = colors_comb(iLocComb_all(iLoc), :);
        %         bar(iData+(iLoc-2)*buffer, CoV_med(iData, iLoc), 'BarWidth', buffer, 'FaceColor', 'w', 'EdgeColor', color)
        %         errorbar(iData+(iLoc-2)*buffer, CoV_med(iData, iLoc), CoV_SEM_neg(iData, iLoc), CoV_SEM_pos(iData, iLoc), 'CapSize', 0, 'color', color)
        
        histogram(CoV_allB{iData, iLoc}, nBins, 'FaceColor', color, 'FaceAlpha', .3, 'EdgeColor', 'none', 'Normalization', 'probability')
        xline(CoV_med(iData, iLoc), 'color', color);
        
    end % iLoc
    
    ylabel('Coefficient of Variation')
%     yticks(linspace(0, .08, 5)), ylim([0, .08])

    x_ticks = round(linspace(x_ticks_all(iData, 1), x_ticks_all(iData, 2), 5), 2);
    xticks(x_ticks), xlim(x_ticks([1, end]))
    
    title(t_all{iData})
    set(findall(gcf, '-property', 'LineWidth'), 'LineWidth',1.5)
    set(findall(gcf, '-property', 'fontsize'), 'fontsize', 18)
    
    saveas(gcf, sprintf('%s/%s/CoV/%s.jpg', nameFigFolder, name_numFilters_Fitting, t_all_short{iData}))
    
end %iData
% xticks(1:nData)
% xticklabels(t_all), xtickangle(45)



