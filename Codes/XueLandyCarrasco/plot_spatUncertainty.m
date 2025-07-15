
subjList = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA', 'CS', 'DT'};
nblocks_allSubj = [200, 240, 210, 240, 220, 220, 215, 220, 205, 210, 205, 205];
nsubj=length(subjList);
SX_RC1_setting

wd_lr = 1;
nLoc5 = 5;

y_itcpt = nan(nsubj, 1); x_itcpt =y_itcpt;
x_uniq_allSubj = nan(nsubj, iLoc);
y_uniq_allSubj = x_uniq_allSubj;
figure('Position', [0 0 1.5e3 1e3]), hold on

for isubj = 1:nsubj
    load(sprintf('Data_OOD/%s%d/%s_behavMeas.mat', subjList{isubj}, nblocks_allSubj(isubj), subjList{isubj}))
    cst_perSess_perLoc = cst_perSess_perLoc.^2;
    x_uniq_allLoc = [];
    y_uniq_allLoc = x_uniq_allLoc;
    
    corr_text = [];
    subplot(3,4, isubj), hold on
    for iLoc = 1:5
        
        x = cst_perSess_perLoc(:, iLoc);
        y = dprime_perSess_perLoc(:, iLoc);
        
        x_uniq=unique(cst_perSess_perLoc(:, iLoc));
        nx=length(x_uniq);
        y_uniq=nan(nx, 2);
        
        for ix=1:nx
            dprime_all = dprime_perSess_perLoc(cst_perSess_perLoc(:, iLoc)==x_uniq(ix), iLoc);
            y_uniq(ix, 1)=mean(dprime_all);
            y_uniq(ix, 2)=std(dprime_all);
        end
        x_uniq_allLoc = [x_uniq; x_uniq_allLoc];
        y_uniq_allLoc = [y_uniq; y_uniq_allLoc];
        errorbar(x_uniq, y_uniq(:, 1), y_uniq(:, 2), 'o', 'capsize', 0, 'color', colors_comb(iLoc, :))
        %         plot(x_uniq, y_ave, 'o')
        
        % linear regression and corr
        lm = polyfit(x_uniq, y_uniq(:, 1), 1);
        %         x_lm2 = linspace(min(x_uniq), max(x_uniq), 2);
        x_lm2 = linspace(-.1, max(x_uniq), 2);
        yfit = polyval(lm, x_lm2);
        %         plot(x_lm2, yfit,'-', 'color',colors_comb(iLoc, :) , 'handlevisibility', 'off', 'linewidth', wd_lr);
        eta2_all = var(polyval(lm, x_uniq))/var(y_uniq(:, 1));
        
        [r, p] = corr(x_uniq, y_uniq(:, 1));
        corr_text = [corr_text, sprintf('\nL%d: r=%.2f, p=%.3f', iLoc, r, p)];
        
        
        x_uniq_allSubj(isubj, iLoc) = mean(x_uniq);
        y_uniq_allSubj(isubj, iLoc) = mean(y_uniq(:,1));
    end % iLoc
    
    
    % linear regression and corr
    lm = polyfit(x_uniq_allLoc, y_uniq_allLoc(:, 1), 1);
    x_lm2 = linspace(min(x_uniq_allLoc), max(x_uniq_allLoc), 2);
    x_lm2 = linspace(-.1, max(x_uniq_allLoc), 2);
    yfit = polyval(lm, x_lm2);
    plot(x_lm2, yfit,'-k' , 'handlevisibility', 'off', 'linewidth', wd_lr*3);
    eta2_all = var(polyval(lm, x_uniq_allLoc)/var(y_uniq_allLoc(:, 1)));
    [r, p] = corr(x_uniq_allLoc, y_uniq_allLoc(:, 1));
    corr_text = [corr_text, sprintf('\nr=%.2f, p=%.3f', r, p)];
    
    %         xlim([-.6,-.1])
    xlim([-.1,.8]), xticks(0:.2:.8)
    ylim([-.5, 2.5]), yticks(0:2)
    
    yline(0, 'color', ones(1,3)/2);
    xline(0, 'color', ones(1,3)/2);
    %
    %     xticks(-.6:.1:-.1)
    %     xticklabels(round(10.^(-.6:.1:-.1),2))
    xintercept = -lm(2)/lm(1);
    
    
    title(sprintf('%s%s', subjList{isubj}, corr_text))
    title(sprintf('%s\ny-intercept=%.2f', subjList{isubj}, yfit(1)))
    title(sprintf('S%d\nx-intercept=%.2f, y-intercept=%.2f', isubj, xintercept, yfit(1)))
    title(sprintf('S%d', isubj))
    y_itcpt(isubj) = yfit(1);
    x_itcpt(isubj) = xintercept(1);
    if isubj==1
        xlabel('Contrast threshold')
        ylabel('dprime')
    end
end

% set(findall(gcf, '-property', 'LineWidth'), 'LineWidth',2)
set(findall(gcf, '-property', 'fontsize'), 'fontsize',20)
%%
figure('Position', [0 0 1.5e3, 400])
subplot(1,2,1), hold on
for iLoc = 1:nLoc5
    for isubj=1:nsubj
        plot(x_uniq_allSubj(isubj, iLoc), y_uniq_allSubj(isubj, iLoc),markers_allSubj{isubj} , 'Color', colors_comb(iLoc, :), 'MarkerSize', 10)
    end
end

lm = polyfit(x_uniq_allSubj(:), y_uniq_allSubj(:), 1);
x_lm2 = linspace(0, max(x_uniq), 2);
yfit = polyval(lm, x_lm2);
plot(x_lm2, yfit,'-k' , 'handlevisibility', 'off', 'linewidth', wd_lr);

subplot(1,2,2), hold on
% [y_itcpt_ave,~, ~, y_itcpt_SEM] = getCI(y_itcpt, 2, 1);
% bar(1, y_itcpt_ave)
% errorbar(1, y_itcpt_ave,y_itcpt_SEM, 'CapSize', 0 )

boxplot([x_itcpt , y_itcpt],'Labels',{'x-intercept','y-intercept'})
[~, px, ~, statsx] = ttest(x_itcpt);
[~, py, ~, statsy] = ttest(y_itcpt);
yline(0, '-k');
title(sprintf('x: t=%.2f, p=%.3f\ny: t=%.2f, p=%.3f', statsx.tstat, px, statsy.tstat, py))

set(findall(gcf, '-property', 'LineWidth'), 'LineWidth',2)
set(findall(gcf, '-property', 'fontsize'), 'fontsize',30)

