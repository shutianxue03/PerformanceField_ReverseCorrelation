% Model Recovery
% estimate internal-external noise ratio by minimizing the nLL calculated
% based on trial-wise responses
% version 2: estimate parameters by search thru the grid

% Observation:
% (1) The simulation and estimation across iterations are the same (weird)
% (2) The simulation and prediction are very similar (good), but still slightly away from the mathematical result from simulation
% (3) The estimated params are still off...

addpath(genpath('Data_PR_OOD'))
addpath(genpath('fxn_model'))
addpath(genpath('fxn_analysis_RC'))

clc
close all

%% simulation SETTING
truth.mu_PRS = 0.64; % mean of the internal variable distribution [PRS]
truth.sigma_PRS = 0.09; % std of the internal variable distribution [PRS], treated as the external noise
truth.mu_ABS = 0.53;
truth.sigma_ABS = 0.08;
params_true = [1.5, 1.2, .6];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ntrials_sim = 1e5; % [simulation] num of trials simulated to generate the data set
truth.nparams_model = 3;
errorComp = 1:8;
fitMode = 2; % (1) fit 8 metrics (2) fit trial-wise responses

%% model fitting SETTING
% namesMetrics_short = {'dprime',  'Criterion', 'pC', 'pHit', 'pFA', 'pA', 'pA [PRS]', 'pA [ABS]'};
% nmetrics = length(namesMetrics_short);
% options = optimoptions('fmincon','MaxIterations',5000,'Display','off');
% options_GS = optimoptions(@fmincon,'Algorithm','sqp');

%% simulation
sim = PR_sim(ntrials_sim, truth, params_true, 0);

%%
n1 = 20; %has to be odd number
n2 = 20;
n3 = 20;
% set all possible values each parameter could take
% a. the true value is  NOT listed in the vector (naturally)
param1_all = [linspace(0, 2, n1)];
param2_all = [linspace(0, 2, n2)];
param3_all = [linspace(0, 2, n3)];

% b. the true value is listed in the vector
param1_all = [linspace(0,2, n1-1)];  param1_all = [param1_all(1:14), params_true(1), param1_all(15:end)];
param2_all = [linspace(0, 2, n2-1)]; param2_all = [param2_all(1:11), params_true(2), param2_all(12:end)];
param3_all = [linspace(0, 2, n3-1)]; param3_all = [param3_all(1:6), params_true(3), param3_all(7:end)];

param_all_all = {param1_all, param2_all, param3_all};

% get the grid value closest to the true value
[~, iparam1_grid] = min(abs(param1_all - params_true(1)));
[~, iparam2_grid] = min(abs(param2_all - params_true(2)));
[~, iparam3_grid] = min(abs(param3_all - params_true(3)));
iparam_grid_all = [iparam1_grid, iparam2_grid, iparam3_grid];

% get nLL for the true value
nLL_true = fxn_getError(fitMode, params_true, sim, truth, errorComp);

% get nLL for the grid value closest to the true value
nLL_grid = fxn_getError(fitMode, [param1_all(iparam1_grid), param2_all(iparam2_grid), param3_all(iparam3_grid)], sim, truth, errorComp);

%% get nLL of each param comb
nLL_3D = fxn_estParams_grid(sim, truth, fitMode, param_all_all, errorComp);

%% get the index of the min
error_min = min(nLL_3D(:));
for ip1 = 1:n1 % do not change to combvec, parfor does not accept that
    for ip2 = 1:n2
        for ip3 = 1:n3
            error_picked = nLL_3D(ip1, ip2, ip3);
            % compare the error picked with the minimum
            if abs(error_picked - error_min)< 1e-5, iparams_est = [ip1, ip2, ip3]; end
        end
    end
end

params_est = [param1_all(iparams_est(1)), param2_all(iparams_est(2)), param3_all(iparams_est(3))];

%% plot
close all

xlabels = {'alpha [PRS]', 'alpha [ABS]', 'thresh'};

figure('Position', [2000 0 1500 400])
for ip = 1:3
    param_all = param_all_all{ip};
    iparam_grid = iparam_grid_all(ip); % the point on grid closest to the true value
    
    subplot(1,3,ip), hold on
    stem(param_all, param_all, 'handlevisibility', 'off')
    xline(params_true(ip), 'r-', 'linewidth', 2); % true value
    stem(param_all(iparam_grid), param_all(iparam_grid), 'r') % the point on grid closest to the true value
    xline(params_est(ip), 'b-', 'linewidth', 2);
    
    if ip == 1, legend({'True param', 'Grid point closest to truth', 'estimated param'}, 'Location', 'south'), end
    xlabel(xlabels{ip})
    xticks([])
    title(sprintf('True param: %.1f\nGrid: [%dth] %.4f\nEstimation: [%dth] %.4f', ...
        params_true(ip), ...
        iparam_grid, param_all(iparam_grid), ...
        iparams_est(ip), params_est(ip)))
end

set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)

sgtitle(sprintf('ntrials for simu: %d\nnLL for the true params: %.2f\nnLL for the grid params: %.2f\nnLL for the current estimation: %.2f\nnLL grid vs. nLL est: %.4f', ...
    ntrials_sim, ...
    nLL_true, ...
    nLL_grid, ...
    error_min, ...
    log(abs(nLL_true-nLL_grid)/(nLL_true+nLL_grid))))

%% 2D plot at each slice

iaxis = [2,1,3; 3,2,1; 3,1,2];
buffer_all = [.02, .002, .002];
iparam = 3; % z is (1) thresh (2) alpha PRS (3) alpha ABS

close all

x_all = param_all_all{iaxis(iparam, 1)}; 
y_all = param_all_all{iaxis(iparam, 2)}; 
z_all = param_all_all{iaxis(iparam, 3)};
buffer = (max(x_all)-min(x_all))/10;
xtruth = params_true(iaxis(iparam, 1));
ytruth = params_true(iaxis(iparam, 2));
ztruth = params_true(iaxis(iparam, 3));
xlabel_ = xlabels{iaxis(iparam, 1)};
ylabel_ = xlabels{iaxis(iparam, 2)};
title_ = xlabels{iaxis(iparam, 3)};

nz = length(z_all);
trajectory = cell(nz, 4);

for it = 1:nz
    switch iparam
        case 1, nLL_2D = squeeze(nLL_3D(:,:,it));
        case 2, nLL_2D = squeeze(nLL_3D(it,:,:));
        case 3, nLL_2D = squeeze(nLL_3D(:,it,:));
    end
    min_nLL = min(nLL_2D(:));
    [min_y, min_x] = find(abs(min_nLL - nLL_2D)<eps);
    
    figure('Position', [1000 500 500 500])
    ylim([min(y_all), max(y_all)])
    xlim([min(x_all), max(x_all)])
    
    % plot 2D
    if abs(z_all(it) - ztruth) < 1e-5, colors = 'w'; else,  colors = 'r'; end
    hold on
    imagesc(x_all, y_all, nLL_2D)
    plot(x_all(min_x), y_all(min_y), [colors, '*']) % local min
    text(x_all(min_x)-buffer, y_all(min_y), sprintf('%.4f', z_all(it)), 'HorizontalAlignment', 'center', 'color', 'r')
    % save local min (to draw a trajectory)
    trajectory(it, :) = {x_all(min_x), y_all(min_y), z_all(it), colors};
    for iit = 1:(it-1)
        plot(trajectory{iit, 1}, trajectory{iit, 2}, [trajectory{iit, 4}, '*'])
        text(trajectory{iit, 1}-buffer, trajectory{iit, 2}, sprintf('%.4f', trajectory{iit, 3}), 'HorizontalAlignment', 'center', 'color', trajectory{iit, 4})
    end
    
    colorbar
%     caxis([0, 1.3e5])
    xline(xtruth, [colors,'-'], 'linewidth', 1.5);
    yline(ytruth, [colors,'-'], 'linewidth', 1.5);
    xlabel(xlabel_)
    ylabel(ylabel_)
    zlabel('nLL')
    title(sprintf('%s = %.3f (%.3f)', title_, z_all(it), ztruth))
    waitforbuttonpress
    %                 pause(.5)
end
