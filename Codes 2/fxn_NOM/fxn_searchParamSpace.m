
%% sim
sim = PR_sim(ntrials_sim, truth, params_true, 0);

%% estimate the param by searching thru the grid
params_est_grid = fxn_estParams_grid(sim, truth, errorComp);

%%
pred = PR_pred(params_est_grid, truth, sim);
metric_pred = pred.metrics;
    
%% 2D plot at each slice

iaxis = [2,1,3; 3,2,1; 3,1,2]; 
buffer_all = [.02, .002, .002];
for iparam = 3 % z is (1) thresh (2) alpha PRS (3) alpha ABS
    
    close all

    xall = p_all{iaxis(iparam, 1)}; yall = p_all{iaxis(iparam, 2)}; zall = p_all{iaxis(iparam, 3)};
    buffer = (max(xall)-min(xall))/10;
    xtruth = params_true_ratio(iaxis(iparam, 1));
    ytruth = params_true_ratio(iaxis(iparam, 2));
    ztruth = params_true_ratio(iaxis(iparam, 3));
    xlabel_ = titles_ratio{iaxis(iparam, 1)};
    ylabel_ = titles_ratio{iaxis(iparam, 2)};
    title_ = titles_ratio{iaxis(iparam, 3)};
    
end
nall = length(zall);

traj = cell(nall, 4);

for it = 1:nall
    switch iparam
        case 1, error2D = squeeze(error_3D(:,:,it));
        case 2, error2D = squeeze(error_3D(it,:,:));
        case 3, error2D = squeeze(error_3D(:,it,:));
    end
    min_error = min(error2D(:));
    [min_y, min_x] = find(min_error == error2D);
      
    figure('Position', [1000 500 500 500])
    ylim([min(yall), max(yall)])
    xlim([min(xall), max(xall)])
    % plot 2D
    if abs(zall(it) - ztruth) < 1e-5, colors = 'w'; else,  colors = 'r'; end
    hold on
    imagesc(xall, yall, error2D)
    plot(xall(min_x), yall(min_y), [colors, '*']) % local min
    text(xall(min_x)-buffer, yall(min_y), sprintf('%.4f', zall(it)), 'HorizontalAlignment', 'center', 'color', 'r')
    % save local min (to draw a trajectory)
    traj(it, :) = {xall(min_x), yall(min_y), zall(it), colors};
    for iit = 1:(it-1)
        plot(traj{iit, 1}, traj{iit, 2}, [traj{iit, 4}, '*'])
        text(traj{iit, 1}-buffer, traj{iit, 2}, sprintf('%.4f', traj{iit, 3}), 'HorizontalAlignment', 'center', 'color', traj{iit, 4})
    end
    
    colorbar
    xline(xtruth, [colors,'-'], 'linewidth', 1.5);
    yline(ytruth, [colors,'-'], 'linewidth', 1.5);
    xlabel(xlabel_)
    ylabel(ylabel_)
    zlabel('nLL')
    title(sprintf('%s = %.3f (%.3f)', title_, zall(it), ztruth))
    waitforbuttonpress
    %                 pause(.5)
end


