

RCplot_setting
if nsubj==1, errType = 1; else, errType = 2; end %1=68%CI, 2=SEM
save('params_RCplot_temp', 'n*', 'names*', 'axis*', 'colors*','errType', 'subjName', 'subplot_locs', 'combInd', 'dprime_theo');

%% plot behav measurements
% group average and SEM
bm_aveSubj = squeeze(mean(data_aveBlk_allSubj, 1)); % where is this defined??
bm_semSubj = squeeze(std(data_aveBlk_allSubj, [], 1))/sqrt(nsubj);
bm_perSess_all = {};
if nsubj==1
    bm_perSess_all = {...
        dprime_perSess_perLoc, ...
        criterion_perSess_perLoc, ...
        RT_perSess_perLoc, ...
        log10(cst_perSess_perLoc), ...
        squeeze(pA_perSess_perLoc(3, :,:)), ...
        squeeze(pA_perSess_perLoc(1, :,:)), ...
        squeeze(pA_perSess_perLoc(2, :,:)), ...
        squeeze(pC_perSess_perLoc(3, :,:)), ...
        squeeze(pC_perSess_perLoc(1, :,:)), ...
        squeeze(pC_perSess_perLoc(2, :,:))};
end

% settings
title_all = {'dprime', 'criterion', 'RT', 'CST Thresh (log scale)', ...
    'Resp consistency (pA)', 'pA [PRS]', 'pA [ABS]', ...
    'Prop. of correctness', 'Hit rate', 'False alarm rate'};
yline_all = {dprime_theo, 0, nan, nan, nan, nan, nan, .7, .5, .5};
ylim_all = {[.3, 1.7], [-1, 1], [.05, .15], log10([.25, .65]), [.5, 1], [.5, 1], [.5, 1], [.5, 1], [.5, 1], [0, 1]};
ticks_scatter_all = {0.4:.7:1.8, -1:1:1,  -3:-1, -.5:.15:-.2, .5:.1:1, .5:.1:1, .5:.1:1, .5:.1:1, .5:.1:1, 0:.5:1};
lim_scatter_all = {[.26, 1.94], [-1.2, 1.2], [-3.2,-.8],  [-.57, -.13],[.48, 1.02],[.48, 1.02], [.48, 1.02], [.48, 1.02], [.48, 1.02], [-.1, 1.1]};

for im = 1:length(title_all)
    
    % settings
    title_ = title_all{im};
    if nsubj==1, bm_perSess = bm_perSess_all{im}; end
    yline_ = yline_all{im};
    ylim_ = ylim_all{im};
    ticks_scatter = ticks_scatter_all{im};
    lim_scatter = lim_scatter_all{im};
    
    % IDVD data
    bm_ave_perSubj = squeeze(data_aveBlk_allSubj(:, im, :)); % cst, dprime, criterion, pA, RT
    bm_std_perSubj = squeeze(data_semBlk_allSubj(:, im, :)) * sqrt(nblocks);
    if nsubj == 1
        bm_ave_perSubj = bm_ave_perSubj';
        bm_std_perSubj = bm_std_perSubj';
    end
    
    % Group average
    bm_aveSubj_ = bm_aveSubj(im, :);
    bm_semSubj_ = bm_semSubj(im, :);
    
    %
    % ==================
    RCplot_behavMeas
    % ==================
end


%% 2. plot 2D kernels and related
if nsubj > 1
    kernels2D = squeeze(median(kernels2D_allSubj, 2));
    intercept2D = squeeze(median(intercept2D_allSubj, 2));
    R2_2D = squeeze(median(R2_2D_allSubj, 2));
    R2_Tjur = squeeze(median(R2_Tjur_allSubj, 2));
    pValues = squeeze(median(pValues_allSubj, 2));
    pCat = squeeze(median(pCat_allSubj, 2));
    seprblity = squeeze(median(sep_allSubj, 2));
else
    kernels2D = squeeze(kernels2D_allSubj);
    intercept2D = squeeze(intercept2D_allSubj);
    R2_2D = squeeze(R2_2D_allSubj);
    R2_Tjur = squeeze(R2_Tjur_allSubj);
    pValues = squeeze(pValues_allSubj);
    pCat = squeeze(pCat_allSubj);
    seprblity = squeeze(sep_allSubj);
end

%% 2D
RCplot_2D(kernels2D, 2, 'Kernels')

%% Intercept
RCplot_2D(intercept2D, 1, 'Intercept')

%% Tjur R2
RCplot_2D(R2_Tjur, 1, 'Tjur R^2')

%% p value
if nB > 1, RCplot_2D(squeeze(pValues(:, :, :, :, :, 2)), 1, 'p values'), end

%% proportion of correct categorization
% RCplot_2D(pCat, 1, 'p(correct categ)')

%% separability
RCplot_separability_5Loc
RCplot_separability_combLoc

%% get IDVD median & CI and Group ave and sem
RCplot_getAVE

%% marginalized kernels per 5 loc
RCplot_MARG5Loc

%% marginalized kernels per combined loc
RCplot_MARGCombLoc

%% estimated params (itype=3, all trials)
RCplot_MARGParamsCombLoc

%% rHIT/rFA as a fxn of binned energy
% the input is still allB
if nsubj == 1
    r = pC_tgt_allBB;
    e = ebin_tgt_allBB;
else
    r = squeeze(nanmedian(pC_tgt_norm_allSubj, 2));
    e = squeeze(nanmedian(ebin_tgt_norm_allSubj, 2));
end
RCplot_energyBins

%% correlation between params and threshold
RCplot_corr

%% compare pA PRS vs. ABS
if nsubj>1, RCplot_comp_pA_prsabs, end

%% analysis on behav data for each observer
if nsubj + nB == 2
    compareLocComb % way 1. by concatenating trials % way 2. by averaging the kernels
    RCplot_SAT
    RCplot_stimPhase
    load(namePatchOutputs)
    RCplot_FFT % input: fft1D_mean_both
    RCplot_CI % input: CI*
    plotFlag_unik = 0; RCplot_linearity % do NOT plot idvd block
    plotFlag_unik = 1; RCplot_linearity  % plot idvd block + grouped by cst
    RCplot_pred_pApC % predict pA and pC based on PTM (Dosher & Lu 2008)
end


fprintf('======== PLOTTING DONE ========')
