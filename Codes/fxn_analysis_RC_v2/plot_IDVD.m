
clc
close all
clear all
warning off
format compact

addpath(genpath('Data_OOD'))
addpath(genpath('fxn_analysis_RC_v2'))
addpath(genpath('VSS2023'))

%%
SX_RC1_setting
subjList = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA', 'CS', 'DT'};
nblocks_allSubj = [200, 240, 210, 240, 220, 220, 215, 220, 205, 210, 205, 205];
nsubj = length(subjList);

%% select trials
% based on two criteria
% 1. pC is within a range
% 2. pA is above a value

% criterion
buffer = .15; ub_pC = .7+buffer; lb_pC = .7-buffer; % pC
lb_pA = .55; % pA has to be > 0.55 (0.5 is the chance level)
iSess_start = 6; % start from the 6th session

% empty containers
bm_perSess_allSubj = cell(1, nsubj);
iSess_delete_pC_allSubj = cell(1, nsubj); iSess_select_pC_allSubj = cell(1, nsubj);
iSess_delete_pA_allSubj = cell(1, nsubj); iSess_select_pA_allSubj = cell(1, nsubj);
iSess_delete_allSubj = cell(1, nsubj); iSess_select_allSubj = cell(1, nsubj);
usability_allSubj = nan(nsubj, 3);

for isubj = 1:nsubj
    subjName = subjList{isubj};
    nblocks = nblocks_allSubj(isubj);
    
    load(sprintf('data_OOD/%s%d/%s_behavMeas.mat', subjName, nblocks, subjName));
    nSess = size(dprime_perSess_perLoc,1);
    nAllTrials = (nSess- iSess_start+1)*100;
    
    bm_perSess_all = {...
        dprime_perSess_perLoc, ...
        criterion_perSess_perLoc, ...
        squeeze(pC3_perSess_perLoc(:, :, 1)), ...
        squeeze(pC3_perSess_perLoc(:, :, 2)), ...
        squeeze(pC3_perSess_perLoc(:, :, 3)), ...
        squeeze(pA3_perSess_perLoc(:, :, 1)), ...
        squeeze(pA3_perSess_perLoc(:, :, 2)), ...
        squeeze(pA3_perSess_perLoc(:, :, 3)), ...
        1./cst_perSess_perLoc, ...
        log(RT_perSess_perLoc)};
    
    bm_perSess_allSubj{isubj}= bm_perSess_all;
    
    %% criterion 1: pC is within a range
    iSess_select_pC = []; iSess_delete_pC = [];
    iLoc_delete_pC = zeros(nSess- iSess_start+1,nLoc5); % which loc does not meet the criterion
    for iSess = iSess_start:nSess
        if (sum(bm_perSess_all{3}(iSess, :) > ub_pC)) || (sum(bm_perSess_all{3}(iSess, :) < lb_pC))
            iSess_delete_pC = [iSess_delete_pC, iSess];
            for iLoc = 1:nLoc5
                if (bm_perSess_all{3}(iSess, iLoc) > ub_pC) || (bm_perSess_all{3}(iSess, iLoc) < lb_pC)
                    iLoc_delete_pC(iSess, iLoc) = 1;
                end
            end
        else
            iSess_select_pC = [iSess_select_pC, iSess];
        end % iSess
    end
    assert(length(iSess_delete_pC) + length(iSess_select_pC) == nSess - iSess_start+1)
    
    %% criterion 2: pA is above a value
    iSess_select_pA = []; iSess_delete_pA = [];
    iLoc_delete_pA = zeros(nSess- iSess_start+1,nLoc5); % which loc does not meet the criterion
    for iSess = iSess_start:nSess
        if sum(bm_perSess_all{6}(iSess, :) < lb_pA)
            iSess_delete_pA = [iSess_delete_pA, iSess];
            for iLoc = 1:nLoc5
                if bm_perSess_all{6}(iSess, iLoc) < lb_pA
                    iLoc_delete_pA(iSess, iLoc) = 1;
                end
            end
        else
            iSess_select_pA = [iSess_select_pA, iSess];
        end % iSess
    end
    assert(length(iSess_delete_pA) + length(iSess_select_pA) == nSess- iSess_start+1)
    
    %%
    iSess_delete = unique([iSess_delete_pC, iSess_delete_pA]);
    iSess_select = 1:nSess; iSess_select(iSess_delete) = []; iSess_select = iSess_select(iSess_select>iSess_start-1); assert(length(iSess_delete) + length(iSess_select) == nSess - iSess_start+1)
    iSess_delete_pC_allSubj{isubj} = iSess_delete_pC;
    iSess_select_pC_allSubj{isubj} = iSess_select_pC;
    iSess_delete_pA_allSubj{isubj} = iSess_delete_pA;
    iSess_select_pA_allSubj{isubj} = iSess_select_pA;
    iSess_delete_allSubj{isubj} = iSess_delete;
    iSess_select_allSubj{isubj} = iSess_select;
    %%
    pSelect = length(iSess_select)/(nSess-iSess_start+1);
    usability_allSubj(isubj, :) = [round(pSelect*100), round(nAllTrials*pSelect), nAllTrials];
    fprintf('%s: %d%% usable, %d/%d trials per Loc pC [%s] %d pA [%s] %d\n', subjName, usability_allSubj(isubj, :), ...
        num2str(sum(iLoc_delete_pC)), sum(iLoc_delete_pC(:)), num2str(sum(iLoc_delete_pA)), sum(iLoc_delete_pA(:)))
end % isubj

data = [round(mean(usability_allSubj)); round(std(usability_allSubj))];
fprintf('\nSummary: nsubj = %d, %d%% (%d%%) usable, %d(%d)/%d(%d) trials per Loc\n', nsubj, data(:))
save('Data_OOD/n12_data_usability', '*select*allSubj*', '*delete*allSubj*', 'usability_allSubj')

%% Fig 1. plot the behavMeasures of each session for all 5 locs
% select subj
SX_RC1_setting
subjList = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA', 'CS', 'DT'};
nblocks_allSubj = [200, 240, 210, 240, 220, 220, 215, 220, 205, 210, 205, 205];
eyeD_all = [1,1,0,1,1,1,1,1,0,0,1,0]; % 1=right eye dominant; 0=left eye dominant

indSubj = [1:6, 8, 9, 12]; % no CS, HA, AS
indSubj = 1:12;
nsubj = length(indSubj); fprintf('nsubj = %d\n', nsubj)
subjList_ = subjList(indSubj);
nblocks_allSubj = nblocks_allSubj(indSubj);
markers_allSubj = markers_allSubj(indSubj);

data = [round(mean(usability_allSubj(indSubj, :))); round(std(usability_allSubj(indSubj, :)))];
fprintf('\nSummary: nsubj = %d, %d%% (%d%%) usable, %d(%d)/%d(%d) trials per Loc\n', length(indSubj), data(:))

%%%%%%%%%%%%%%%%
namesMetrics_plus2 = [namesMetrics, 'CS', 'RT'];

ticks_scatter_all = {0:.5:3, -1:.5:1, 0:.2:1, 0:.2:1, 0:.2:1, 0:.2:1, 0:.2:1, 0:.2:1, 0:4, -5:1:-1};
lim_scatter_all = {[0, 3], [-1, 1],  [.5, 1], [0, 1], [0, 1], [.5, 1], [0, 1], [0, 1], [1, 4], [-5,-1]};
yline_all = {dprime_theo, 0, threshPerf, .5, .5, .5, .5, .5, nan, nan};

for isubj = 1:nsubj
    subjName = subjList_{isubj};
    nblocks = nblocks_allSubj(isubj);
    bm_perSess_all = bm_perSess_allSubj{isubj};
    iSess_select_pC = iSess_select_pC_allSubj{isubj};
    iSess_select_pA = iSess_select_pA_allSubj{isubj};
    iSess_select = iSess_select_allSubj{isubj};
    iSess_delete_pC = iSess_delete_pC_allSubj{isubj};
    iSess_delete_pA = iSess_delete_pA_allSubj{isubj};
    iSess_delete = iSess_delete_allSubj{isubj};
    nSess = nblocks/nLoc5;
    
    for im = 1:nmetrics+2
        bm_perSess = bm_perSess_all{im};
        figure('Position', [0 500 2000 300])
        
        %%%%%%%%%%%%%%
        % change
        %%%%%%%%%%%%%%
        subplot(1,6,1:4), hold on, grid on
        for iLoc = 1:nLoc5
            if iLoc== 4, linestyle = '--'; else, linestyle = '-'; end
            plot(bm_perSess(:,iLoc), linestyle, 'color',colors_comb(iLoc, :), 'LineWidth', .5),
        end
        % extra line
        if ~isnan(yline_all{im}), yline(yline_all{im}, '-', 'color', ones(1,3)*.5, 'linewidth', 2); end
        %
        for iSess = 1:nSess
            if sum(iSess == iSess_delete_pA), xline(iSess, ':k', 'linewidth', 2); end
            if sum(iSess == iSess_delete_pC), xline(iSess, '--k', 'linewidth', 2); end
        end
        ylim(lim_scatter_all{im})
        yticks(ticks_scatter_all{im})
        if im == 10, yticklabels(round(exp(ticks_scatter_all{im})*100+500)), end
        
        xlabel('Session #')
        ylabel(namesMetrics_plus2{im})
        xtickInd = round(quantile(1:nSess, linspace(0,1,5)));
        if nSess > 1, xlim([0, nSess + 2])
            if nSess>2, xticks(xtickInd), xticklabels(xtickInd), else, xticks([1,2]), end
        else, xlim([.5,2]), xticks(1)
        end
        
        % references (the criterion of exclusion)
        if im == 3
            patch([[0, nSess + 2], flip([0, nSess + 2])], [[lb_pC, lb_pC], [ub_pC, ub_pC]], 'r', 'FaceAlpha', .1, 'linestyle', 'none')
        elseif im==6
            yline(lb_pA, 'r-', 'linewidth', 2);
        end
        
        switch im
            case 3, yticks(.5:.1:1)
            case 6, yticks(.5:.1:1)
            case 9, yticks(1:4)
        end
        
        %%%%%%%%%%%%%%
        %  bar 1 (raw data)
        %%%%%%%%%%%%%
        subplot(1,6,5), hold on, grid on
        fxn_bars(im, nLoc5, bm_perSess, colors_comb, yline_all, namesLocComb, lim_scatter_all, ticks_scatter_all)
        
        %%%%%%%%%%%%%%
        %  bar 2 (selected data)
        %%%%%%%%%%%%%
        subplot(1,6,6), hold on, grid on
        fxn_bars(im, nLoc5, bm_perSess(iSess_select, :), colors_comb, yline_all, namesLocComb, lim_scatter_all, ticks_scatter_all)
        
        %
        set(findall(gcf, '-property', 'FontSize'), 'FontSize',18)
        pSelect = length(iSess_select)/nSess;
        sgtitle(sprintf('%s: %d%% usable, %d/%d trials per Loc\n', subjName, usability_allSubj(indSubj(isubj), :)), 'FontSize',22)
        
        folderName = sprintf('VSS2023/fig/dataUsability/%s/', namesMetrics_plus2{im});
        folderDir = dir(folderName);
        if isempty(folderDir), mkdir(folderName), end
        
        saveas(gcf, sprintf('%s%s.jpg', folderName, subjName))
        
    end % im
end % isubj

%% Fig 2. plot cst and pC in one figure
% RULE: if pC drops, cst should also drop
nCol=7;
lim_pC = [.55, .85];
lim_cst = [.2, .7];

for isubj = 1:nsubj
    figure('Position', [0 0 2000 800])
    subjName = subjList{isubj};
    nblocks = nblocks_allSubj(isubj);
    bm_perSess_all = bm_perSess_allSubj{isubj};
    nSess = length(bm_perSess_all{1});
    for iLoc = 1:nLoc5
        subplot(nLoc5, nCol, (1+nCol*(iLoc-1)):(nCol*iLoc-3)), hold on
        % pC
        yyaxis left
        data = bm_perSess_all{3}(:, iLoc);
        data_select = data(5:end); % all sessions after the first 5
        data_smoothed = smoothdata(data, 'gaussian', 5);
        plot(1:nSess, data_smoothed, '-')
        ylabel('pC')
        yline(.7, 'k-');
        ylim(lim_pC)
        xlim([0, 49])
        xticks(5:5:45)
        %         errorbar(44.5, mean(data), std(data), 'bo', 'CapSize', 0) % average of all sessions
        %         errorbar(45, mean(data_select), std(data_select), 'bo', 'CapSize', 0) % average of all sessions after the first 5
        % CST
        yyaxis right
        data_cst = 1./bm_perSess_all{9}(:, iLoc);
        data_cst_select = data_cst(end-10:end); % last 10 sessions
        plot(1:nSess, data_cst, '-')
        ylim(lim_cst)
        ylabel('CST')
        %         errorbar(45.5, mean(data_cst), std(data_cst), 'ro', 'CapSize', 0) % average of all sessions
        %         errorbar(46, mean(data_cst_select), std(data_cst_select), 'ro', 'CapSize', 0) % average of last 10 sessions
        
        title(namesLocComb{iLoc})
    end % iLoc
    
    % compare pC/cst across loc
    for im = [3,9]
        if im == 3, data = bm_perSess_all{im}; data_select = data(5:end, :); name = 'pC';% pC: all sessions after the first 5
        else, data = 1./bm_perSess_all{im}; data_select = data(end-10:end, :); name = 'CST'; % cst: last 10 sessions
        end
        for ii = 1:5
            if im == 3, subplot(nLoc5, nCol, nCol*ii-1), ylim(lim_pC), yline(.7, 'k-'); % pC
            else, subplot(nLoc5, nCol, nCol*ii), ylim(lim_cst) % cst
            end
            hold on
            switch ii
                case 1 % all 5 loc
                    for iLoc = 1:5
                        errorbar(iLoc, mean(data(:, iLoc)), std(data(:, iLoc)), 'o', 'color', colors_comb(iLoc, :), 'CapSize', 0)
                        errorbar(iLoc+6, mean(data_select(:, iLoc)), std(data_select(:, iLoc)), 'o', 'color', colors_comb(iLoc, :), 'CapSize', 0)
                    end
                    title(sprintf('[%s] 5 loc', name))
                    xlim([0, 12])
                    xline(6, '-', 'color', ones(1,3)*.7);
                    xticks([])
                case 2 % F vs. P
                    fxn_plot(data, data_select, 1, 2:5, colors_comb(1, :), colors_comb(8, :))
                    title(sprintf('[%s] Fov vs. Peri', name))
                    xlim([0,6])
                case 3 % L vs. R
                    fxn_plot(data, data_select, 2, 4, colors_comb(2, :), colors_comb(4, :))
                    title(sprintf('[%s] Left vs. Right', name))
                    xlim([0,6])
                case 4 % HM vs. VM
                    fxn_plot(data, data_select, [2,4], [5,3], colors_comb(6, :), colors_comb(7, :))
                    title(sprintf('[%s] HM vs. VM', name))
                    xlim([0,6])
                case 5 % LVM vs. UVM
                    fxn_plot(data, data_select, 5,3, colors_comb(5, :), colors_comb(3, :))
                    title(sprintf('[%s] LVM vs. UVM', name))
                    xlim([0,6])
                    xticks([1.5, 4.5])
                    xticklabels({'all sess', 'from 6th sess'}), xtickangle(45)
            end % switch ii
        end % ii
    end % im
    
    set(findall(gcf, '-property', 'fontsize'), 'fontsize', 15)
    set(findall(gcf, '-property', 'LineWidth'), 'LineWidth',1.5)
    sgtitle(subjName)
    
    saveas(gcf, sprintf('VSS2023/fig/dataUsability/_pC_and_CST/%s.jpg', subjName))
    
end % isubj

close all

%% helper fxn - fxn_bars
function fxn_bars(im, nLoc5, bm_perSess, colors_comb, yline_all, namesLocComb, lim_scatter_all, ticks_scatter_all)
for iLoc = 1:nLoc5
    bar(iLoc, mean(bm_perSess(:, iLoc)), 'FaceColor', 'w', 'EdgeColor', colors_comb(iLoc, :))
    errorbar(iLoc, median(bm_perSess(:, iLoc)), std(bm_perSess(:, iLoc)), '.', 'color', colors_comb(iLoc, :), 'CapSize', 0)
end

% extra line
if ~isnan(yline_all{im}), yline(yline_all{im}, '-', 'color', ones(1,3)*.5, 'linewidth', 2); end

xticks(1:nLoc5)
xticklabels(namesLocComb(1:nLoc5))
ylim(lim_scatter_all{im})
yticks(ticks_scatter_all{im})
yticklabels([])
end


%% helper fxn - fxn_plot
function fxn_plot(data, data_select, iLoc_all1, iLoc_all2, color1, color2)
getRatio = @(a,b) (a-b)./(a+b);

data1 = data(:,iLoc_all1); data1 = data1(:);
data2 = data(:,iLoc_all2); data2 = data2(:);
yyaxis left
errorbar(1, mean(data1),  std(data1), 'o', 'color', color1, 'CapSize', 0)
errorbar(2, mean(data2),  std(data2), 'o', 'color', color2, 'CapSize', 0)
yyaxis right
bar(1.5, getRatio(mean(data1), mean(data2)), 'EdgeColor', 'k', 'Facecolor', 'w', 'barwidth', .2)
yticks(-.1:.1:.1)
ylim([-.1, .1])

data1_s = data_select(:,iLoc_all1); data1_s = data1_s(:);
data2_s = data_select(:,iLoc_all2); data2_s = data2_s(:);
yyaxis left
errorbar(4, mean(data1_s),  std(data1_s), 'o', 'color', color1, 'CapSize', 0)
errorbar(5, mean(data2_s),  std(data2_s), 'o', 'color', color2, 'CapSize', 0)
yyaxis right
bar(4.5, getRatio(mean(data1_s), mean(data2_s)), 'EdgeColor', 'k', 'Facecolor', 'w', 'barwidth', .2)
yticks(-.1:.1:.1)
ylim([-.1, .1])
ylabel('Ratio')
xline(3, '-', 'color', ones(1,3)*.7);
xticks([])
end
