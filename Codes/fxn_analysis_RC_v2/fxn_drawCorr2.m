function  p = fxn_drawCorr2(x_allSubj, y_allSubj, xticks_, yticks_, ANOVA_indLoc, iLocComb_all, flag_plotCI, flag_plotINSET,flag_plotConnect,  markers_allSubj, colors_comb)
%% figure setting
wd_border = 6;
sz_ticks = 50; % 30
sz_marker = 40;


nsubj = size(x_allSubj, 1);
if ndims(x_allSubj)==2, nLoc = size(x_allSubj, 2);
else, nLoc = size(x_allSubj, 3);
end


%% zero-mean
x_zeroM = x_allSubj-mean(x_allSubj, 1);
y_zeroM = y_allSubj - mean(y_allSubj, 1);

%%
figure('Position', [0 200 1200 1000]); hold on, box on
% 
% %% connect idvd data
% if flag_plotConnect
%     if nLoc==2
%         for isubj = 1:nsubj
%             % if the estP/tunC of HM/LVM > VM/UVM, solid line; otherwise, dashed line
%             if y_zeroM(isubj, 1) > y_zeroM(isubj, 2), linestyle = '-'; else, linestyle = ':'; end
%             plot(x_zeroM(isubj, :), y_zeroM(isubj, :), ['k.', linestyle], 'linewidth', wd_border/1.5)
%         end
%     end
% end

%% idvd markers with CI
for iiLoc = 1:nLoc
    for isubj = 1:nsubj
%         if flag_plotCI
%             errorbar(x_zeroM(isubj, iiLoc), y_zeroM(isubj, iiLoc), ...
%                 y_neg(isubj, iiLoc), y_pos(isubj, iiLoc), x_neg(isubj, iiLoc), x_pos(isubj, iiLoc), ...
%                 'color', colors_comb(iLocComb_all(iiLoc), :), 'linewidth', wd_border, 'capsize', 0)
%         end
        plot(x_zeroM(isubj, iiLoc), y_zeroM(isubj, iiLoc),  markers_allSubj{isubj}, 'markerfacecolor', 'w', 'markeredgecolor', colors_comb(iLocComb_all(iiLoc), :), 'markersize', sz_marker, 'linewidth', wd_border)
        %         plot(x_zeroM(isubj, iLoc), y_zeroM(isubj, iLoc),  'o', 'markerfacecolor', 'w', 'markeredgecolor', colors_comb(iLocComb_all(iLoc), :), 'markersize', sz_marker, 'linewidth', wd_border)
    end % isubj
end % iLoc

%% get partial corr
[r_partial, p_partial] = partialcorr(x_allSubj(:), y_allSubj(:), ANOVA_indLoc);
% get eta2 (effect size)

%% corr for each loc
text_corr_perL = [];
for iiLoc = 1:nLoc
    [r, p] = corr(x_allSubj(:, iiLoc), y_allSubj(:, iiLoc));
    text_corr_perL = [text_corr_perL, sprintf('L%d: r=%.2f, p=%.3f\n', iLocComb_all(iiLoc), r, p)];
end

%% linear regression
lm = polyfit(x_zeroM(:), y_zeroM(:), 1);
x_lm2 = linspace(min(x_zeroM(:)), max(x_zeroM(:)), 2);
yfit = polyval(lm, x_lm2);
plot(x_lm2, yfit,'-', 'color', ones(1,3)*.4, 'handlevisibility', 'off', 'linewidth', wd_border * 1.5);
eta2_all = var(polyval(lm, x_zeroM(:)))/var(y_zeroM(:));

%%

% xlim_ = xticks_([1, end]);
% ylim_ = yticks_([1, end]);
% xlim(xlim_), xticks(xticks_) % temporaily hard-coded, since x is usually CS
% ylim(ylim_), yticks(yticks_)

axis square
title(sprintf('Partial r = %.2f (p = %.3f) eta^2=%.2f\n%s', round(r_partial,2), round(p_partial,3), round(eta2_all, 2), text_corr_perL))

ax = gca;
ax.XAxis.FontSize = sz_ticks;
ax.YAxis.FontSize = sz_ticks;
ax.LineWidth = wd_border;

%% inset (raw data)
% if ip==3, flag_plotINSET=1;end
if flag_plotINSET
    axes('Position',[.53 .6 .3 .3])
    
    box on, hold on
    for iiLoc = 1:nLoc
        cc = colors_comb(iLocComb_all(iiLoc), :);
        for isubj = 1:nsubj
            plot(x_allSubj(isubj, iiLoc), y_allSubj(isubj, iiLoc), markers_allSubj{isubj}, ...
                'MarkerSize', sz_marker/1.5, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', cc, 'linewidth', wd_border/1.5)
        end
    end
    
    %     xlim(xlim_inset), xticks(xticks_inset)
    %     ylim(ylim_inset), yticks(yticks_inset)
    
    axis square
    ax = gca;
    ax.XAxis.FontSize = sz_ticks/1.5;
    ax.YAxis.FontSize = sz_ticks/1.5;
    ax.LineWidth = wd_border/1.5;
end % if flag_plotINSET

