clc
varNames = {'Loc', 'Channel'};

wd_border = 8;
sz_ticks = 45; % 30

shift = [0,0];

iLoc_all_all = {[6,7], [5,3]};

for ii = 1:2
    iLoc_all = iLoc_all_all{ii};
    if ii == 1, peak = [.3, .2; .3, .28];
    else, peak = [.2, .1; .2, .18];
    end
    
colors_comb_ = colors_comb(iLoc_all, :);
figure('Position', [0 0 7e2 1e3])
    for ifeature = 1:2
        
        xaxis = axis_tuning{ifeature};
        
        %% plot kernels & tuning fxn
        subplot(2, 1, ifeature)
        hold on, box on
        % extra lines
        yline(0, 'handlevisibility', 'off', 'linewidth', wd_border/2, 'color', [.7, .7, .7]);
        xline(ifeature-1, 'handlevisibility', 'off', 'linewidth', wd_border/2, 'color', [.7, .7, .7]);
        for iLoc = 1:nLoc
            color = colors_comb_(iLoc, :);
            if ifeature == 1
                if iLoc == 1, margPred_ave = predSFkernel(xaxis, 1, [peak(ifeature, iLoc) 17, 0], 0);
                else, margPred_ave = predSFkernel(xaxis, 1, [peak(ifeature, iLoc), 17, 0], 0);
                end
            else
                if iLoc == 1, margPred_ave = predSFkernel(2.^xaxis, 2, [1.6, peak(ifeature, iLoc), .3, 0], 0);
                else, margPred_ave = predSFkernel(2.^xaxis, 2, [1.6, peak(ifeature, iLoc), .3, 0], 0);
                end
            end
            
            if nLoc == 2
                if nsubj>1
                    plot(xaxis, margPred_ave', '-', 'color', color, 'linewidth', wd_border)
                else
                    plot(xaxis, margPred_ave, '-', 'color', color, 'linewidth', wd_border)
                end
            end
        end % end of iiLoc
        
        % xaxis
        xlim(axisLim{ifeature})
        xticks(axisTicks_tuning{ifeature})
        xticklabels(axisTL_tuning{ifeature})
        yticks([0])
        ylim([-.05, .35])
        
        %%
        ax = gca;
        ax.XAxis.FontSize = sz_ticks;
        ax.YAxis(1).FontSize = sz_ticks;
        ax.YAxis(1).Color= 'k';
        ax.LineWidth = wd_border/1.5;
        
        
    end % end of ifeature
    saveas(gcf, sprintf('L%d%d.jpg', iLoc_all))
end

