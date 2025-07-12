
% model simulation
% We explored a predefined set of parameter values for each model variations (iModelB) in order to find the
% upper and lower bound of the parameter space
% iModelB=1: lapse rate, Nmul, SDadd, thresh
% iModelB=2: lapse rate, SDadd, thresh
% iModelB=3: lapse rate, Nmul, thresh
% iModelB=4: lapse rate, thresh

addpath(genpath('fxn_NOM'))
clc, close all
load('NOM_data')
iModelB = 1;

% [.001, .2, .2, 1] gives reasonable prediction

ni_lapse = 1; i_lapse_all = .02; % lapse rate
ni_Nmul = 30; i_Nmul_all = linspace(.00001, .5, ni_Nmul); % Nmul
ni_SDadd = 30; i_SDadd_all = linspace(.00001, .5, ni_SDadd); % SDadd
ni_crit = 30; i_crit_all = linspace(.5, 2.5, ni_crit); % threshold, as the range of IVs is around 0-1.2
switch iModelB
    case 1, nLL = nan(ni_lapse, ni_Nmul, ni_SDadd, ni_crit); ind = combvec(i_lapse_all, i_Nmul_all, i_SDadd_all, i_crit_all);
    case 2, nLL = nan(ni_lapse, ni_SDadd, ni_crit); ind = combvec(i_lapse_all, i_SDadd_all, i_crit_all);
    case 3, nLL = nan(ni_lapse, ni_Nmul, ni_crit); ind = combvec(i_lapse_all, i_Nmul_all, i_crit_all);
    case 4, nLL = nan(ni_lapse, ni_crit); ind = combvec(i_lapse_all, i_crit_all);
end
nInd = size(ind, 2);

%%
mu_IV_PRS = median(data.IV{1});
mu_IV_ABS = median(data.IV{2});
SD_IV_PRS = std(data.IV{2});
SD_IV_ABS = std(data.IV{2});

figure
hold on
histogram(data.IV{1}, 'facecolor', 'r')
histogram(data.IV{2}, 'facecolor', 'b')
legend({'Signal-present', 'Signal-absent'})
xline(mu_IV_PRS, 'r', 'linewidth', 2);
xline(mu_IV_ABS, 'b', 'linewidth', 2);
xlabel('Internal variable')

%%
close all
for iInd = 1:nInd
    if iInd<=10, flag_plot=1; else, flag_plot=0; end
    flag_plot=0;
    params_est = ind(:, iInd);
    switch iModelB
        case 1
            i_lapse = find(params_est(1) == i_lapse_all);
            i_Nmul = find(params_est(2) == i_Nmul_all);
            i_SDadd = find(params_est(3) == i_SDadd_all);
            i_crit = find(params_est(4) == i_crit_all);
            nLL(i_lapse, i_Nmul, i_SDadd, i_crit) = fxn_getError_v5(iModelB, params_est, data, flag_plot);
            
        case 2
            i_lapse = find(params_est(1) == i_lapse_all);
            i_SDadd = find(params_est(2) == i_SDadd_all);
            i_crit = find(params_est(3) == i_crit_all);
            nLL(i_lapse, i_SDadd, i_crit) = fxn_getError_v5(iModelB, params_est, data, flag_plot);
            
        case 3
            i_lapse = find(params_est(1) == i_lapse_all);
            i_Nmul = find(params_est(2) == i_Nmul_all);
            i_crit = find(params_est(3) == i_crit_all);
            nLL(i_lapse, i_Nmul, i_crit) = fxn_getError_v5(iModelB, params_est, data, flag_plot);
    end
end

nLL = squeeze(nLL);
fprintf('\n\nDONE\n\n')

%%
nLL_min = min(nLL(:));
[min_row, min_col] = find(nLL_min == nLL); % min_x is

figure, hold on
switch iModelB
    case 1
        for ii=1:length(i_Nmul_all)
            nLL_slice = squeeze(nLL(ii, :, :));
            nLL_min = min(nLL_slice(:));
            [min_row, min_col] = find(nLL_min == nLL_slice); % min_x is
            imagesc(i_SDadd_all, i_crit_all, nLL_slice), colorbar
            plot(i_SDadd_all(min_col), i_crit_all(min_row), 'r*')
            xlabel('SDadd'), xlim(i_SDadd_all([1,end]))
            ylabel('Thresh'), ylim(i_crit_all([1,end]))
            title(sprintf('Nmul=%.2f (%d/%d)', i_Nmul_all(ii), ii, length(i_Nmul_all)))
            pause(.1)
        end
        
    case 2
        imagesc(i_SDadd_all, i_crit_all, nLL), colorbar
        plot(i_SDadd_all(min_col), i_crit_all(min_row), 'r*')
        xlabel('SDadd'), xlim(i_SDadd_all([1,end]))
        ylabel('Thresh'), ylim(i_crit_all([1,end]))
        
    case 3
        imagesc(i_Nmul_all, i_crit_all, nLL), colorbar
        plot(i_Nmul_all(min_col), i_crit_all(min_row), 'r*')
        xlabel('Nmul'), xlim(i_Nmul_all([1,end]))
        ylabel('Thresh'), ylim(i_crit_all([1,end]))
end


