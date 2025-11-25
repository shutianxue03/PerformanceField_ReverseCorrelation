% plotAll_part1.m

% Last updated by Shutian Xue on 07/22/2025
%
% Description:
%   This script generates summary and visualization plots for reverse correlation analysis.
%   It loads compiled and bootstrapped kernel data, produces figures for 2D kernels, marginalized tuning curves, tuning characteristics, and inter-parameter correlations.
%
% Usage:
%   Run after SX_analysis4_Boot.m to visualize and summarize RC results for all observers

% Helper functions:
%   RCplot_2Dkernel (function)
%   plot5_marg_grant (script)
%   plot5b_marg_params_idvd (script)
%   plot6_params_grant (script)

close all, clc

%--------------%
SX_RC1_setting
%--------------%

% Create directory to save figures

nameFolderCompile = nameFolder_compile;
% nameFolder_Figures_RC = 'XueCarrasco_JN/fig/RC';

nameFolder_Figures_RC = sprintf('%s/RC_%d%d%s_ORI%dSF%d', nameFolder_Figures, nORI, nSF, text_cut,  ifamily_perF);
if ~exist(nameFolder_Figures_RC, 'dir'), mkdir(nameFolder_Figures_RC); end

% settings
paramMode=2;
flag_standEnergy = 1;
flag_mirrorMapping = 1;
iType = 2;
nLoc = 2;
if nLoc == 2, nameFileLoc = sprintf('L%d%d', iLocComb_all(1), iLocComb_all(2)); else, nameFileLoc = sprintf('all%d', nLoc); end
% name_numFilters_Fitting = sprintf('%d_%d_ORI%d_SF%d', nORI, nSF, ifamily_perF);
% if flag_mirrorMapping == 1
%     name_numFilters_Fitting = sprintf('%d_%d_m_ORI%d_SF%d', nORI, nSF, flag_mirrorMapping, ifamily_perF);
% end
    

%% 1. 2D kernels
clc, close all
e2D_med_allSubj = getCI(kernels2D_allSubj, 1, 2); % med kernel of each subj (to conduct t-test on FP diff)

    nameFolder_kernels2D = sprintf('%s/kernels2D', nameFolder_Figures_RC);
    if isempty(dir(nameFolder_kernels2D)), mkdir(nameFolder_kernels2D), end

% if flag_standEnergy
%     kernelLim_allType{1} = {[-.06, .12], [-.06, .12], [-.04, .08]}; % PRS
kernelLim_allType{2} = {[-.06, .14], [-.06, .14], [-.05, .05]}; % ABS; Loc1, Loc2, diff
%     kernelLim_allType{3} = {[-.05, .1], [-.05, .1], [-.05, .1]}; % BOTH
% else, kernelLim_allType = {[-5, 8], [-5, 8], [-3.5, 5]};
% end

sz_title = 60;
sz_all = 55;
wd_all = 3;

alpha = .01;

for iType = 2
    kernelLim = kernelLim_allType{iType};
    pCriterion = .001;
    [e2D_ave, ~, ~, e2D_sd] = getCI(e2D_med_allSubj(:, :, iType, :, :), 2, 1);

    %%% NEW BOOTS %%%
    %     e2D_aveSubj_allB = getCI(kernels2D_allSubj, 2, 1); % nsubj x nBoot x nLoc x nTypes x nORI x nSF
    %     [e2D_ave, e2D_aveSubj_lb, e2D_aveSubj_ub] = getCI(e2D_aveSubj_allB(:, :, itype, :, :), 1, 1);

    % t-test per channel
    for iORI = 1:nORI
        for iSF = 1:nSF
            [~, p_pos1] = ttest(squeeze(e2D_med_allSubj(:, 1, iType, iORI, iSF)), 0, 'tail', 'right');
            [~, p_pos2] = ttest(squeeze(e2D_med_allSubj(:, 2, iType, iORI, iSF)), 0, 'tail', 'right');

            p_pos1_all(iORI, iSF) = p_pos1*nORI*nSF;
            p_pos2_all(iORI, iSF) = p_pos2*nORI*nSF;
        end
    end
    [y_interp, x_interp] = meshgrid(axis_tuning{2}, -90:90);

    e2D_ave1_interp = nan(nSF, 181);
    e2D_ave2_interp = e2D_ave1_interp;
    for iSF =1:nSF
        e2D_ave1_interp(iSF, :) =  interp1(axis_tuning{1}, squeeze(e2D_ave(1, :, iSF)), -90:90);
        e2D_ave2_interp(iSF, :) =  interp1(axis_tuning{1}, squeeze(e2D_ave(2, :, iSF)), -90:90);
    end
    e2D_diff_interp = e2D_ave1_interp - e2D_ave2_interp;

    fprintf('[%s]\n', namesType{iType})
    for iiLoc=1:2
        ss = squeeze(e2D_ave(iiLoc, :, :));
        fprintf('%s: MIN is %.3f, MAX is %.3f\n', namesLocComb{iLocComb_all(iiLoc)}, min(ss(:)), max(ss(:)))
    end

    %%%%%%
    % Loc 1
    %%%%%%
    figure('Position', [0 0 1e3 1e3])
    %     RCplot_2Dkernel(e2D_ave1_interp, kernelLim{1}, p_pos1_all'<=alpha)
    RCplot_2Dkernel(e2D_ave1_interp, kernelLim{1})
    %     title(sprintf('[%s] %s', namesType{itype}, namesLocComb{iLocComb_all(1)}), 'FontSize', sz_title)
    title(namesLocComb{iLocComb_all(1)}, 'FontSize', sz_title)
    saveas(gcf, sprintf('%s/n%d_L%d_%s.jpg', nameFolder_kernels2D, nsubj, iLocComb_all(1), namesType{iType}))

    %%%%%%
    % Loc 2
    %%%%%%
    figure('Position', [4e2 0 1e3 1e3])
    %     RCplot_2Dkernel(e2D_ave2_interp, kernelLim{2}, p_pos2_all'<=alpha)
    RCplot_2Dkernel(e2D_ave2_interp, kernelLim{2})
    %     title(sprintf('%s', namesType{itype}, namesLocComb{iLocComb_all(2)}), 'FontSize', sz_title)
    title(namesLocComb{iLocComb_all(2)}, 'FontSize', sz_title)
    saveas(gcf, sprintf('%s/n%d_L%d_%s.jpg', nameFolder_kernels2D, nsubj, iLocComb_all(2), namesType{iType}))

    %%%%%%
    % Diff
    %%%%%%
    figure('Position', [8e2 0 1e3 1e3])
    % diff = squeeze(e2D_med_allSubj(:, 1, :, :) - e2D_med_allSubj(:, 2, :, :));
    %     diff_ave = squeeze(e2D_ave(1, :, :) - e2D_ave(2, :, :));
    fprintf('DIFF: MIN is %.2f, MAX is %.2f\n', min(e2D_diff_interp(:)), max(e2D_diff_interp(:)))
    h_pos = nan(nORI, nSF); h_neg = h_pos;
    for iORI = 1:nORI
        for iSF = 1:nSF
            %             h_pos(iORI, iSF) = ttest(e2D_med_allSubj(:, 1, itype, iORI, iSF), e2D_med_allSubj(:, 2, itype, iORI, iSF),  'alpha', .05/nORI/nSF, 'tail', 'right');
            h_pos(iORI, iSF) = ttest(e2D_med_allSubj(:, 1, iType, iORI, iSF), e2D_med_allSubj(:, 2, iType, iORI, iSF),  'alpha', .05/nORI/nSF);
            %             h_neg(iORI, iSF) = ttest(e2D_med_allSubj(:, 1, itype, iORI, iSF), e2D_med_allSubj(:, 2, itype, iORI, iSF),  'alpha', .05/nORI/nSF, 'tail', 'left');
        end
    end
    RCplot_2Dkernel(e2D_diff_interp, kernelLim{3}, h_pos')
    %     title(sprintf('%s minus %s', namesType{itype}, namesLocComb{iLocComb_all(1)}, namesLocComb{iLocComb_all(2)}), 'FontSize', sz_title)
    title(sprintf('%s minus %s', namesLocComb{iLocComb_all(1)}, namesLocComb{iLocComb_all(2)}), 'FontSize', sz_title)
    saveas(gcf, sprintf('%s/n%d_L%d%d_diff_%s.jpg', nameFolder_kernels2D, nsubj, iLocComb_all(1), iLocComb_all(2), namesType{iType}))
end
% close all

%% 1b. 2D kernels for each idvd
e2D_med_allSubj = getCI(kernels2D_allSubj, 1, 2); % med kernel of each subj (to conduct t-test on FP diff)

    nameFolder_kernel2D_idvd = sprintf('%s/%s/idvd/kernels2D/L%d%d/', nameFolder_Figures_RC, name_numFilters_Fitting, iLocComb_all);
    if isempty(dir(nameFolder_kernel2D_idvd)), mkdir(nameFolder_kernel2D_idvd), end

[y_interp, x_interp] = meshgrid(axis_tuning{2}, -90:90);
for isubj = 1:nsubj
    %     figure('Position', [3e3 0 1.2e3 1.2e3])
    figure('Position', [0 0 1.2e3 4e2])
    for iType = 2
        e2D_med = getCI(e2D_med_allSubj(isubj, :, iType, :,: ), 2, 1);

        % interpolate because the cell size are constant even for fine-sampling
        e2D_med1_interp = nan(nSF, 181);
        e2D_med2_interp = e2D_med1_interp;
        for iSF =1:nSF
            e2D_med1_interp(iSF, :) =  interp1(axis_tuning{1}, squeeze(e2D_med(1, :, iSF)), -90:90);
            e2D_med2_interp(iSF, :) =  interp1(axis_tuning{1}, squeeze(e2D_med(2, :, iSF)), -90:90);
        end
        e2D_diff_interp = e2D_med1_interp - e2D_med2_interp;

        %%%%%%
        % Loc 1
        %%%%%%
        %         subplot(3,3,1+(itype-1)*3), hold on
        subplot(1,3,1), hold on
        RCplot_2Dkernel(e2D_med1_interp, kernelLim{1})
        set(findall(gcf, '-property', 'FontSize'), 'FontSize',18)
        title(sprintf('[%s] %s', namesType{iType}, namesLocComb{iLocComb_all(1)}), 'FontSize',20)

        %%%%%%
        % Loc 2
        %%%%%%
        %         subplot(3,3,2+(itype-1)*3), hold on
        subplot(1,3,2), hold on
        RCplot_2Dkernel(e2D_med2_interp, kernelLim{2})
        set(findall(gcf, '-property', 'FontSize'), 'FontSize',18)
        title(sprintf('[%s] %s', namesType{iType}, namesLocComb{iLocComb_all(2)}), 'FontSize',20)

        %%%%%%
        % Diff
        %         %%%%%%
        %         subplot(3,3,3+(itype-1)*3), hold on
        subplot(1,3,3), hold on
        RCplot_2Dkernel(e2D_diff_interp, kernelLim{3})
        title(sprintf('[%s] %s minus %s', namesType{iType}, namesLocComb{iLocComb_all(1)}, namesLocComb{iLocComb_all(2)}), 'FontSize',20)

    end % itype
    set(findall(gcf, '-property', 'FontSize'), 'FontSize',18)
    sgtitle(subjList{isubj})

    % save
    saveas(gcf, sprintf('%s%s.jpg', folderName, subjList{isubj}))
end
close all

%% still 2D mapping, but all subj on the same figure, each figure is for one loc
iType = 2;
for iiLoc=1:2
    figure('Position', [0, 0, 2e3, 1.8e3])
    for isubj = 1:nsubj

        e2D = squeeze(e2D_med_allSubj(isubj, iiLoc, iType, :,:));

        % interpolate because the cell size are constant even for fine-sampling
        e2D_med_interp = nan(nSF, 181);
        for iSF =1:nSF
            e2D_med_interp(iSF, :) =  interp1(axis_tuning{1}, squeeze(e2D(:, iSF)), -90:90);
        end

        subplot(3,5,isubj), hold on
        RCplot_2Dkernel(e2D_med_interp, kernelLim{1})
        title(subjList{isubj})

    end % isubj
    set(findall(gcf, '-property', 'FontSize'), 'FontSize',12)
    sgtitle(sprintf('[%s] %s', namesType{iType}, namesLocComb{iLocComb_all(iiLoc)}), 'FontSize',20)

    % save
    folderName = sprintf('%s/%s/idvd/kernels2D/L%d%d/', nameFolder_Figures_RC, name_numFilters_Fitting, iLocComb_all);
    folderDir = dir(folderName); if isempty(folderDir), mkdir(folderName), end
    saveas(gcf, sprintf('%sL%d_allSubj.jpg', folderName, iLocComb_all(iiLoc)))
end % iLoc
% close all

%% 5. marginalized kernels per combined loc (compare two loci)
% clc
flag_plotSFPeak = 0; % 1=plot SF peaks;
flag_plotIDVD = 0;
xTLs = namesLocComb(iLocComb_all);
title_ = sprintf('L%d%d', iLocComb_all);

for iType=2%:3
    %     plot5_marg
    plot5_marg_grant
    %     plot5_marg_pred
end

%% 5b. [IDVD tuning] for each loc pair, plot the median & 68% CI of kernels & diff
% if nsubj==15
for iType = 2
    for flagPlotAllSubjInOnePlot = 1%[1,0]
        % 1=plot ONLY the kernels&tuning (% diff) of all subj in one plot
        % 0=plot for each subj: kernels, tuning, diff, params
        for ifeature = 1:2
            plot5b_marg_params_idvd
        end
    end
end % itype
% end % if nsubj==12
% close all

%% 6. params (all trials)
close all, clc
flag_plotIDVD = 0;
flag_plotDiff = 1;

plot6_params_grant

%% reorganize two peaks based on dominance
if ifamily_perF(2)==12
    margTuningC_SF_allSubj_ = margTuningC_SF_allSubj;

    turns = zeros(nsubj, nB, 2);
    iType=2;
    for isubj=1:12
        for iiLoc = 1:2
            for iB=1:nB
                % variation 1: first 3 belong to the left peak, last 3 belong to the right peak, with 4th being the distance to 2 cpd
                %                 margTuningC_SF_allSubj_(isubj, iB, iiLoc, itype, 4) = abs(log2(margTuningC_SF_allSubj(isubj, iB, iiLoc, itype, 4))-1);

                % variation 2: LAST 3 belongs to the dominant peak (with the higher peakAmp), FIRST 3 belong to non-dom peak
                %                 if margTuningC_SF_allSubj(isubj, iB, iiLoc, itype, 2) >= margTuningC_SF_allSubj(isubj, iB, iiLoc, itype, 5)
                %                     margTuningC_SF_allSubj_(isubj, iB, iiLoc, itype, 4) = margTuningC_SF_allSubj(isubj, iB, iiLoc, itype, 1);
                %                     margTuningC_SF_allSubj_(isubj, iB, iiLoc, itype, 5) = margTuningC_SF_allSubj(isubj, iB, iiLoc, itype, 2);
                %                     margTuningC_SF_allSubj_(isubj, iB, iiLoc, itype, 6) = margTuningC_SF_allSubj(isubj, iB, iiLoc, itype, 3);
                %
                %                     margTuningC_SF_allSubj_(isubj, iB, iiLoc, itype, 1) = margTuningC_SF_allSubj(isubj, iB, iiLoc, itype, 4);
                %                     margTuningC_SF_allSubj_(isubj, iB, iiLoc, itype, 2) = margTuningC_SF_allSubj(isubj, iB, iiLoc, itype, 5);
                %                     margTuningC_SF_allSubj_(isubj, iB, iiLoc, itype, 3) = margTuningC_SF_allSubj(isubj, iB, iiLoc, itype, 6);
                %                     turns(isubj, iB, iiLoc)=1;
                %                 end

                % variation 3: based on var 2: with 1st and 4th being the distance to 2 cpd
                if margTuningC_SF_allSubj(isubj, iB, iiLoc, iType, 2) >= margTuningC_SF_allSubj(isubj, iB, iiLoc, iType, 5)
                    margTuningC_SF_allSubj_(isubj, iB, iiLoc, iType, 4) = abs(log2(margTuningC_SF_allSubj(isubj, iB, iiLoc, iType, 1))-1);
                    margTuningC_SF_allSubj_(isubj, iB, iiLoc, iType, 5) = margTuningC_SF_allSubj(isubj, iB, iiLoc, iType, 2);
                    margTuningC_SF_allSubj_(isubj, iB, iiLoc, iType, 6) = margTuningC_SF_allSubj(isubj, iB, iiLoc, iType, 3);

                    margTuningC_SF_allSubj_(isubj, iB, iiLoc, iType, 1) = abs(log2(margTuningC_SF_allSubj(isubj, iB, iiLoc, iType, 4))-1);
                    margTuningC_SF_allSubj_(isubj, iB, iiLoc, iType, 2) = margTuningC_SF_allSubj(isubj, iB, iiLoc, iType, 5);
                    margTuningC_SF_allSubj_(isubj, iB, iiLoc, iType, 3) = margTuningC_SF_allSubj(isubj, iB, iiLoc, iType, 6);
                else
                    margTuningC_SF_allSubj_(isubj, iB, iiLoc, iType, 1) = abs(log2(margTuningC_SF_allSubj(isubj, iB, iiLoc, iType, 1))-1);
                    margTuningC_SF_allSubj_(isubj, iB, iiLoc, iType, 4) = abs(log2(margTuningC_SF_allSubj(isubj, iB, iiLoc, iType, 4))-1);
                end
            end% iB
        end%iiLoc
    end % isubj
end


%% 6a. when SF Family is 2/3, plot SF peak of each idvd
% xTLs = namesLocComb_;
% title_ = sprintf('L%d%d', iLocComb_all);

plot5_marg_grant

%% 6b. when SF Family is 12, which SF peak the system gives more weight to
wd = 5;
wd_bar = 2;
sz_marker_ave = 30;
sz_marker_idvd = 12;
sz_ticks = 20;

if ifamily_perF(2) == 12
    peakAmp = getCI(margTuningC_SF_allSubj(:, :, :, iType, [2,5]), 1, 2); % nsubj x nLoc x [peakSF1, peakSF2]
    ind_Loc = repmat(1:nLoc, nsubj, 1, 2);
    ind_whichPeak = cat(3, ones(nsubj,nLoc), ones(nsubj,nLoc)*2);

    tt = print_nANOVA({'Loc', 'whichPeak'}, peakAmp(:), {ind_Loc(:), ind_whichPeak(:)}, nsubj);
    fprintf('L%d%d:\n%s\n', iLocComb_all, tt)
    [peakAmp_ave, ~, ~, peakAmp_sem] = getCI(peakAmp, 2, 1);

    figure('Position', [0 200 800 300])
    for iLoc = 1:2
        subplot(1,2,iLoc),hold on
        for ipeak = 1:2
            errorbar(ipeak, peakAmp_ave(iLoc, ipeak), peakAmp_sem(iLoc, ipeak), '.', 'color', colors_comb_(iLoc, :), 'CapSize', 0, 'linewidth', wd)
            plot(ipeak, peakAmp_ave(iLoc, ipeak), 'o', 'MarkerEdgeColor', colors_comb_(iLoc, :), 'MarkerSize', sz_marker_ave, 'linewidth', wd)
        end
        buffer = .2;
        for isubj = 1:nsubj
            plot([1+buffer, 2-buffer], squeeze(peakAmp(isubj, iLoc, :)),  [markers_allSubj{isubj}, '-'], ...
                'color',ones(1,3)*.7, 'markerfacecolor', 'w', 'markeredgecolor', ones(1,3)*.7, ...
                'markersize', sz_marker_idvd, 'linewidth', wd_bar)
        end
        [~, p, ~, stats] = ttest(squeeze(peakAmp(:, iLoc, 1)), squeeze(peakAmp(:, iLoc, 2)));
        title(sprintf('t(%d)=%.2f, p=%.3f\n%d/%d', stats.df, stats.tstat, p*2, sum(peakAmp(:, iLoc, 1)>peakAmp(:, iLoc, 2)), nsubj))
        xticks(1:2)
        xticklabels({'Left peak', 'Right peak'})
        xlim([.5, 2.5])
        ylabel('SF peak amplitude')
        yticks(0:.05:.15)
        ylim([0, .15])

        ax = gca;
        ax.XAxis.FontSize = sz_ticks;
        ax.YAxis.FontSize = sz_ticks;
        ax.LineWidth = wd;
    end % iLoc

    % save
    folderName = sprintf('%s/%s/params/%s/%s/', nameFolder_Figures_RC, name_numFilters_Fitting, nameFileLoc, namesType{iType});
    folderDir = dir(folderName); if isempty(folderDir), mkdir(folderName), end
    saveas(gcf, sprintf('%sn%d_comp_two_peaks.jpg', folderName, nsubj))
end


%% 7. assess inter-parameter correlations
% plot7_comp_params_vs_tuningC

% 4 loc
tunC_ORI = margTuningC_ORI_4allSubj;
tunC_SF = margTuningC_SF_4allSubj;
nameFileLoc = 'Lall4';
colors_comb_ = colors_comb(2:5, :);
namesLocComb_ = namesLocComb(2:5);

% 3 loc
tunC_ORI = margTuningC_ORI_3allSubj;
tunC_SF = margTuningC_SF_3allSubj;
nameFileLoc = 'Lall3';
colors_comb_ = colors_comb([6, 5, 3], :);
namesLocComb_ = namesLocComb([6,5,3]);

for ifeature = 1:2
    plot7b_corr_among_tunC
end

% plot7a_comp_estP_tunC
%

