clc
close all

%% loading the compiled data of 5 single locations
% data of all subj (n=12, including CS) are saved in all5, no need to save
% another file for n=11
flag_block200 = 0;
nORI = 29;
nSF = 29;
nB = 1e3; %[nB_n10, nB_nn] = fxn_simplifyNumber(nB);
ifamily_perF = [8,12]; % ORI: Gaussian; SF: log parabola with truncation
flag_PatchMode=2;
flag_standEnergy=1;
text_m = '_m';
if flag_PatchMode == 1 % energy derives from target patch
    nameFileLoc5 = sprintf('Data_compile/T/n12_B%d_Lall4_%d_%d%s_ORI%dSF%d.mat', ...
        nB, nORI, nSF, text_m, ifamily_perF(1), ifamily_perF(2));
    nameEnergySource = 'T';
else % energy derives from noise patch
    nameFileLoc5 = sprintf('Data_compile/N%d/n12_B%d_Lall4_%d_%d%s_ORI%dSF%d.mat', ...
        flag_standEnergy, nB, nORI, nSF, text_m, ifamily_perF(1), ifamily_perF(2));
    nameEnergySource = sprintf('N%d', flag_standEnergy);
end
fprintf('Loading...')
load(nameFileLoc5)
fprintf(' DONE\n')

%%
SX_RC1_setting
clc
subjList = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA', 'CS', 'DT'};
nblocks_allSubj = [200, 240, 210, 240, 220, 220, 215, 220, 205, 210, 205, 205];
eyeD_all = [1,1,0,1,1,1,1,1,0,0,1,0]; % 1=right eye dominant; 0=left eye dominant

indSubj = [1:10, 12]; % no CS
% indSubj = [1:3, 5:10, 12]; % no CS, LS
% indSubj = [1:9, 12]; % no CS, HA
% indSubj = [1:6, 8: 10, 12]; % no CS, AS
% indSubj = 1:12;

nSess_s = 10; % only look at cs from the last 10 sessions
nsubj = length(indSubj); fprintf('nsubj = %d\n', nsubj)
subjList = subjList(indSubj);
nblocks_allSubj = nblocks_allSubj(indSubj);
markers_allSubj = markers_allSubj(indSubj);
eyeD_all = eyeD_all(indSubj);

%% decide location list
indLoc = [2,4,5,3]; %2:nLoc5; % 1:nLoc5 % if plot all 5 locations; 2:nLoc5 if exclude fovea
nLoc = length(indLoc);

%% folder/file names
if flag_block200, nameFigFolder = 'VSS2023/fig_200/RC'; else, nameFigFolder = 'VSS2023/fig/RC'; end
if nLoc == 2, nameFileLoc = sprintf('L%d%d', iLocComb_all(1), iLocComb_all(2)); else, nameFileLoc = sprintf('all%d', nLoc); end

%%
colors_comb = [
    0,0,0; ...,     % center; black
    0, .75, 0; ..., % Left: light green
    1, 0, 0,; ...,   % upper: red
    0, .35, 0; ..., % right: dark green
    0, 0, 1; ...,    % lower: blue
    0, .5, 0; ...,   % HM: green
    .5, 0, 1; ...,    % VM: purple
    .5, .5, .5];      % peri: darker grey

itype = 2;
nLoc = length(indLoc);
namesLocComb_ = namesLocComb(indLoc);
colors_comb_ = colors_comb(indLoc, :);
metrics_allSubj_ = metrics5_allSubj(indSubj, :, indLoc, :);
% RT_allSubj_ = RT5_allSubj(indSubj, :, indLoc);
pYES_tgt_allSubj_ = pYES_tgt5_allSubj(indSubj,:, indLoc, :, :);
ebin_tgt_allSubj_ = ebin_tgt5_allSubj(indSubj,:, indLoc, :, :);
sep_allSubj_ = sep5_allSubj(indSubj, :, indLoc, :);
% marg & pred & R2
margORI_allSubj_ = margORI5_allSubj(indSubj, :, indLoc, :, :);
margPredORI_allSubj_ = margPredORI5_allSubj(indSubj, :, indLoc, :, :);
margR2ORI_allSubj_ = margR2ORI5_allSubj(indSubj, :, indLoc, :);
margSF_allSubj_ = margSF5_allSubj(indSubj, :, indLoc, :, :);
margPredSF_allSubj_ = margPredSF5_allSubj(indSubj, :, indLoc, :, :);
margR2SF_allSubj_ = margR2SF5_allSubj(indSubj, :, indLoc, :);
% estimated params
estP_ORI_allSubj = margParamsORI5_allSubj(indSubj, :, indLoc, :, :);
estP_SF_allSubj = margParamsSF5_allSubj(indSubj, :, indLoc, :, :);
% tuning characteristics
tuningC_ORI_allSubj = margParams2ORI5_allSubj(indSubj, :, indLoc, :, :);
tuningC_SF_allSubj = margParams2SF5_allSubj(indSubj, :, indLoc, :, :);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% recreate cs_allSubj_
% since contrast is not stable until the last several sessions
nameFileCS_selected = sprintf('Data_compile/%s/n%d_B%d_CS_all5.mat', nameEnergySource, nsubj, nB);
% no need to save another file for all4!!
dirCS_selected = dir(nameFileCS_selected);
if isempty(dirCS_selected)
    cst_ave_allSubj = nan(nsubj, nB, 5);
    for isubj = 1:nsubj
        subjName = subjList{isubj};
        
        load(sprintf('Data_OOD/%s%d/%s_behavMeas', subjName, nblocks_allSubj(isubj), subjName))
        cst_perLoc = cst_perSess_perLoc(end-nSess_s+1:end, :);
        % fovea
        cst_ave1 = mean(cst_perLoc(:, 1));
        % left & right
        cst_ave2 = mean(cst_perLoc(:, 2));
        cst_ave4 = mean(cst_perLoc(:, 4));
        % upper and lower
        cst_ave5 = mean(cst_perLoc(:, 5));
        cst_ave3 = mean(cst_perLoc(:, 3));
        
        cst_ave_allSubj(isubj, :, :) = repmat([cst_ave1, cst_ave2, cst_ave3, cst_ave4, cst_ave5], nB, 1);
    end
    
    %%%%%%%%%%%%%
    cs_allSubj_ = 1./cst_ave_allSubj;
    save(nameFileCS_selected, 'cs_allSubj_')
    fprintf('CS file DONE\n')
else, load(nameFileCS_selected), fprintf('CS file Loaded\n')
end
cs_allSubj_ = cs_allSubj_(:, :, indLoc);

%% correlation with CMF
load('VSS2023/CMF_n5.mat') % CMF_n5_area, CMF_n5_subjList
cs_med = getCI(cs_allSubj_, 1, 2);
nLoc_CMF = 4;
markers_CMF = {'o', 's', 'd', '^','v'}; % for each subj
if nsubj == 12
    isubj_CMF = [2, 5, 9, 11, 12]; % SP, RE, FH, CS, DT
end
cs = cs_med(isubj_CMF, 2:5);
nsubj_CMF = length(isubj_CMF); assert(nsubj_CMF == length(CMF_n5_subjList))

figure, hold on
for iisubj_CMF = 1:nsubj_CMF
    for iiLoc_CMF = 1:nLoc_CMF
        plot(cs(iisubj_CMF, iiLoc_CMF), CMF_n5_area(iisubj_CMF, iiLoc_CMF), ...
            markers_CMF{iisubj_CMF}, 'color', colors_comb(iiLoc_CMF+1, :), 'markersize', 10)
    end
end
xlabel('CS')
ylabel('surface area (mm^2)')

[r, p] = corr(cs(:), CMF_n5_area(:));
title(sprintf('RAW DATA\nr=%.3f, p=%.3f', r, p))
set(findall(gcf, '-property', 'LineWidth'), 'LineWidth',2)
set(findall(gcf, '-property', 'fontsize'), 'fontsize',25)

%%
ANOVA_indLoc = repmat(1:nLoc_CMF, nsubj_CMF, 1);
ANOVA_indSubj = repmat((1:nsubj_CMF)', 1, nLoc_CMF);

for dim_zeromean = 1:2 % 1=reveal observer effect; 2=reveal loc effect
    namesZeroMeanDIM = {'reveal observer effect', 'reveal loc effect'};
    
    cs_zeromean = cs - mean(cs, dim_zeromean);
    cmf_zeromean = CMF_n5_area - mean(CMF_n5_area, dim_zeromean);
    figure, hold on
    
    
    for iisubj_CMF = 1:nsubj_CMF
        [~, indOrder] = sort(cs_zeromean(iisubj_CMF, :));
        plot(cs_zeromean(iisubj_CMF, indOrder), cmf_zeromean(iisubj_CMF, indOrder), '-', 'color', ones(1,3)*.5)
        for iiLoc_CMF = 1:nLoc_CMF
            
            plot(cs_zeromean(iisubj_CMF, iiLoc_CMF), cmf_zeromean(iisubj_CMF, iiLoc_CMF), ...
                markers_CMF{iisubj_CMF}, 'markerfacecolor', 'w', 'markeredgecolor', colors_comb(iiLoc_CMF+1, :), 'markersize', 10)
        end
    end
    xlabel('CS')
    ylabel('surface area (mm^2)')
    
    [r, p] = corr(cs_zeromean(:), cmf_zeromean(:));
    if dim_zeromean == 1
        [r, p] = partialcorr([cs(:), CMF_n5_area(:)], ANOVA_indLoc(:)); r = r(2,1); p = p(2,1); title_ = sprintf('Control for loc: Partial r =%.3f, p=%.3f\n', r, p);
    else
        [r, p] = partialcorr([cs(:), CMF_n5_area(:)], ANOVA_indSubj(:)); r = r(2,1); p = p(2,1); title_ = sprintf('Control for subj: Partial r=%.3f, p=%.3f\n', r, p);
    end
    
    
    title(sprintf('%s\nr=%.3f, p=%.3f\n%s', namesZeroMeanDIM{dim_zeromean}, r, p, title_))
    set(findall(gcf, '-property', 'LineWidth'), 'LineWidth', 1.5)
    set(findall(gcf, '-property', 'fontsize'), 'fontsize',20)
end

%%
CI_ratio = .68;

%% 1. plot behav measurements
close all
namesMetrics_plus2 = [namesMetrics, 'CS', 'RT'];
yline_all = {dprime_theo, 0, threshPerf, nan, nan, nan, nan, nan, nan, nan};
ticks_scatter_all = {0:1:2, -1:1:1, .5:.25:1, .5:.25:1, 0:.5:1, .5:.25:1, .5:.25:1, .5:.25:1, 1.5:1:3.5, -4:2:0};
lim_scatter_all = {[0, 2], [-1, 1],  [.5, 1], [.5, 1], [0, 1], [.5, 1], [.5, 1], [.5, 1], [1.5, 3.5], [-4,0]};

flag_plotScatter = 0;
flag_plotIDVD = 0;
% iLocComb_all_ = iLocComb_all;
for im = 1:nmetrics + 1
    plot1_behavMeas
    
end

%% 2. pYES vs. binned energy
plot2_ebins

%% 3. separability
% 4b. get new separability with smaller ORI/SF range
cut_ORI = 4:26;
cut_SF = 6:24;
itype = 2;
axisTicks_tuning = {-60:30:60, [0.3571, 0.6429, 1, 1.3571, 1.6429]}; % ticks (SF is on log scale)
axisTL_tuning = {axisTicks_tuning{1}, round(2.^axisTicks_tuning{2},2)}; % label (SF is on linear scale)
axisLim = {[-99, 99], [-.1, 2.1]};

sep = nan(nsubj, 2);
for isubj = 1:nsubj
    for iiLoc = 1:4
        e2D_cut = squeeze(e2D_med_allSubj(isubj, iiLoc, itype, cut_ORI, cut_SF));
        e2D_recon = (mean(e2D_cut, 1)') * (mean(e2D_cut, 2)');
        e2D_recon = e2D_recon';
        sep(isubj, iiLoc) = corr2(e2D_cut, e2D_recon);
%         figure('Position', [0 0 785 440])
%         subplot(3,5, [2,3,7,8]), imagesc(axis_tuning{2}(cut_SF), axis_tuning{1}(cut_ORI), e2D_cut), axis square, xticks(axisTicks_tuning{2}), xticklabels(round(2.^axisTicks_tuning{2}, 2)), yticks(axisTicks_tuning{1}), yline(0, 'r-'); xline(1, 'r-');
%         subplot(3,5, [4,5,9,10]), imagesc(axis_tuning{2}(cut_SF), axis_tuning{1}(cut_ORI), e2D_recon), axis square, xticks(axisTicks_tuning{2}), xticklabels(round(2.^axisTicks_tuning{2}, 2)), yticks(axisTicks_tuning{1}), yline(0, 'r-'); xline(1, 'r-');
%         subplot(3,5, [1,6]), plot( mean(e2D_cut, 2), axis_tuning{1}(cut_ORI)), set(gca,'YDir','reverse'); yline(0, 'r-');
%         subplot(3,5, [12, 13]), plot(axis_tuning{2}(cut_SF), mean(e2D_cut, 1)), xticks(axisTicks_tuning{2}), xticklabels(round(2.^axisTicks_tuning{2}, 2)), xline(1, 'r-');
%         set(findall(gcf, '-property', 'LineWidth'), 'LineWidth',2)
%         set(findall(gcf, '-property', 'fontsize'), 'fontsize',15)
%         sgtitle(sprintf('%s-L%d (%d%%)', subjList{isubj}, iiLoc, round(sep(isubj, iiLoc)*100)))
    end
end

for itype = 2, plot4_sep, end

%% 5. marginalized kernels per combined loc
title_ = sprintf('all%d', nLoc);
for itype=2, plot5_marg, end

%% 6. compare estP/tunC between/across locations
close all
for itype = 2%:3
flag_plotStats = 1; % 1=plot stats in the title/text
for paramMode=1:2
    if paramMode == 1
        % estimated parameters
        params_allSubj = {estP_ORI_allSubj, estP_SF_allSubj}; title_ = 'estP';
        namesParams = {'peak SF','gain','sigma','baseline','truncation'};
    else
        % tuning characteristics
        params_allSubj = {tuningC_ORI_allSubj, tuningC_SF_allSubj}; title_ = 'tunC';
        namesParams = {'peak SF','peak','width','baseline','truncation'};
    end
    xTLs = namesLoc2D(indLoc);
    plot6_params
end
end % itype

%% 6b. check the corr between corresponding estP/tunC between ORI and SF
plot6b_corr_params_ORI_vs_SF

%% 7. compare estimated params vs. tuning characteristics
plot7a_comp_estP_tunC
plot7b_corr_among_params

%% 8. correlation between params and threshold
for paramMode = 1:2
    switch paramMode
        case 1
            % estimated parameters
            params_allSubj = {estP_ORI_allSubj, estP_SF_allSubj}; title_ = 'estP';
            namesParams = {'peak SF','gain','sigma','baseline','truncation'};
        case 2
            % tuning characteristics
            params_allSubj = {tuningC_ORI_allSubj, tuningC_SF_allSubj}; title_ = 'tunC';
            namesParams = {'peak SF','peak','width','baseline','truncation'};
        case 3
            params_allSubj = {params_med_allSubj, params_allSubj}; title_ = 'NOM';
            namesParams = {'inducedN','constantN','threshold'};
    end
    for itype = 2:3
        for ifeature = 1:2
            plot8_corr
            % adjust im_corr in plot8 to decide which metrics to present
            close all
        end
    end
end

