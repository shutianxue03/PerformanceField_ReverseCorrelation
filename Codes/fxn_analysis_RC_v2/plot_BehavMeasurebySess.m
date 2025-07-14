% Last updated on 07/14/2025 by Shutian Xue
% This scripts  plots the behavioral measures (CST, pA, pC) across sessions for each subject.
% In order to check the stability of the observers' performance across sessions.

iLoc_all = 1:5; % 1=left, 2=right, 3=up, 4=down, 5=center
nLoc = length(iLoc_all);
nMetrics = 3; % 1=CST, 2=pA, 3=pC

nameFolder_Fig_BehavMeasurebySess = sprintf('%s/BehavMeasurebySess', nameFolder_Figures);
if ~exist(nameFolder_Fig_BehavMeasurebySess, 'dir'), mkdir(nameFolder_Fig_BehavMeasurebySess), end

% Preallocate matrices to store measurements
meas_sd = nan(nsubj, nMetrics, nLoc);
meas_min = nan(nsubj, nMetrics, nLoc);
meas_max = nan(nsubj, nMetrics, nLoc);
meas_ave = nan(nsubj, nMetrics, nLoc);

% Loop through each subject
for isubj = 1:nsubj
    
    subjName = subjList{isubj};
    nblock = nblocks_allSubj(isubj);
    
    % Load behavioral measurements for each subject
    load(sprintf('%s/Data_OOD_%d%d/%s%d/%s_behavMeas.mat', nameFolder_Data, nORI, nSF, subjName, nblock, subjName), ...
        'cst_perSess_perLoc', 'pA3_perSess_perLoc', 'pC3_perSess_perLoc')
    
    figure('Position', [0 0 1e3 1e3])
    
    % Loop through each metric (CST, pA, pC)
    for iMetric = 1:nMetrics
        switch iMetric
            case 1
                m= cst_perSess_perLoc*100; 
                y_label = 'Gabor CST (%)'; y_ticks = linspace(30,80, 5); y_ticklabels = y_ticks;
            case 2
                m=squeeze(pA3_perSess_perLoc(:, :, 1))*100; 
                y_label = 'pA (%)'; y_ticks = 40:15:100; y_ticklabels = y_ticks;
            case 3
                m=squeeze(pC3_perSess_perLoc(:, :, 1))*100; 
                y_label = 'pC (%)'; y_ticks = 50:10:90; y_ticklabels = y_ticks;
        end
        subplot(nMetrics,1,iMetric), hold on, grid on
        
        % Loop through each location and plot the measurements
        for iiLoc = 1:nLoc
            plot(m(:, iLoc_all(iiLoc)), 'color',colors_comb(iLoc_all(iiLoc), :)),
        end
        
        % Add reference lines for pA and pC
        if iMetric==2, yline(60, 'k'); yline(80, 'k'); end
        if iMetric==3, yline(70, 'k'); end
        
        ylabel(y_label)
        yticks(y_ticks)
        yticklabels(y_ticklabels)
        ylim(y_ticks([1,end]))
        
        % Calculate and store measurements for each location
        for iiLoc=1:nLoc
            iLoc = iLoc_all(iiLoc);
            meas_sd(isubj, iMetric, iiLoc) = std(m(:, iLoc));
            meas_min(isubj, iMetric, iiLoc) = min(m(:, iLoc));
            meas_max(isubj, iMetric, iiLoc) = max(m(:, iLoc));
            meas_ave(isubj, iMetric, iiLoc) = mean(m(:, iLoc));
        end
        
    end % end of iMetric
    
    xlabel('Session #'),
    xlim([.5, nSess + .5])
    xticks(1:2:nSess)
    %     legend(namesLoc2D, 'location', 'best')
    
    sgtitle(sprintf('%s%d', subjName, nblock))
    
    set(findall(gcf, '-property', 'LineWidth'), 'LineWidth',1.5)
    set(findall(gcf, '-property', 'fontsize'), 'fontsize',15)
    
    % Save figure
    saveas(gcf, sprintf('%s/%s%d.png', nameFolder_Fig_BehavMeasurebySess, subjName, nblock))
end % isubj

%% Conduct ANOVA for each metric across subjects
clc

for iMetric = 3%:nm
    %     for iiLoc = 1:nLoc
    for iStat = 1:4
        switch iStat
            case 1, m = meas_sd; n='SD';
            case 2, m = meas_min; n='MIN';
            case 3, m = meas_max; n='MAX';
            case 4, m = meas_ave; n='AVE';
        end
        
        a=squeeze(m(:, iMetric, :));
        indLoc = repmat(iLoc_all, nsubj, 1);
        text_ANOVA = print_nANOVA({'Loc'}, a(:), {indLoc(:)}, nsubj);
        
        [ave, ~, ~, SEM] = getCI(a(:), 2, 1);
        fprintf('%s: %.1f +- %.1f\n%s\n', n, ave, SEM, text_ANOVA)
    end
end
% CST: [29.7, 52.1%]