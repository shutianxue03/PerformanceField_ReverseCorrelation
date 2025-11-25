

%%%%%%%%%%%%%%%% PART 2 %%%%%%%%%%%%%%%%
% focus on the correlation analysis

close all, clc

addpath(genpath('XueCarrasco_JN/code'))
addpath(genpath('Data_compile/RC'))
addpath(genpath('fxn_analysis_RC_v2'))

%--------------%
SX_RC1_setting
% %--------------%

% fileName_NOM = 'ORI37SF37';
fileName_NOM = 'ORI29SF29';
fprintf('\nnsubj=%d, %s\n', nsubj, fileName_NOM)

% directory
nameFolderCompile_BEHAV = 'Data_compile/BEHAV';
nameFolderCompile_RC = 'Data_compile/RC';
nameFolderCompile_NOM = 'Data_compile/NOM';
nameFigFolder = 'XueCarrasco_JN/fig/RC';

% settings
iLocComb_all = [1,8]; nLoc=3;
colors_comb_ = colors_comb(iLocComb_all, :);
ni = 1e3; % number of interactions (for NOM measures/estimateds)
templateType = 3; % template type: 1=raw, 2=reconstructed, 3=mirrored, 4=m&r
nModelsA = 5; % 1=core, 2-4=random models
iModelA = 1; % 1=only look at core model
paramMode = 2; % 1=estimated params; 2=tuning characteristics
flag_standEnergy = 1;
flag_mirrorMapping = 1;
itype = 2; % 1=PRS, 2=ABS, 3=BOTH

% setting of tuning functions
clc
nORI = 29; nSF = 29; fileName = '2929_m';
name_numFilters_Fitting = sprintf('%d_%d_ORI%d_SF%d', nORI, nSF, ifamily_perF);
if flag_mirrorMapping == 1, name_numFilters_Fitting = sprintf('%d_%d_m_ORI%d_SF%d', nORI, nSF, ifamily_perF); end

nTuningC_ORI = length(namesTunC_unit_perF{ifamily_perF(1), 2}); % preferred ori, peak amp, bandwidth, bottom, // pred_half
nTuningC_SF = length(namesTunC_unit_perF{ifamily_perF(2), 2}); % peakSF, peak amp., bandwidth_full (octave), baseline // bandwidth_L, bandwidth_R, bottom, pred_half, full bandwidth (linear scale)
paramInd_perF{1} = nTuningC_ORI;
paramInd_perF{2} = nTuningC_SF;

%% plot performance
% nSess_s = nan; % select behav of the last 10 sessions (XFC2024); if nan, meaning no select
% fxn_compilePerf_combinedLoc % depends on the number of locations
% fprintf('DONE\n\n')
flag_plotIDVD = 1;
nameFileBehav = 'Lall3';
% nameFileBehav = 'L18';
plotAll_Behav

% fxn_sessSTD % plot how pC and CS change across blocks

%% COMPILE HM, LVM AND UVM (must run thus to plot 3 loc's correlation!!)
clc
if flag_mirrorMapping, text_cut = '_m'; else, text_cut = ''; end
nameFolderCompile = 'Data_compile/RC';
nameFileLoc3 = sprintf('%s/%d%d%s/n%d_B%d_Lall3_%d_%d%s_ORI%dSF%d.mat', ...
    nameFolderCompile, nORI, nSF, text_cut, nsubj, nB, nORI, nSF, text_cut, ifamily_perF);
dirLoc3 = dir(nameFileLoc3);
nLoc3 = 3;

clc, fprintf('\n\nCompiling RC data 3 loc (abs pref ori=%d)...\n', flag_ABSprefORI)

if isempty(dir(nameFileLoc3))
    % empty folders
    similarityAfterMirroring_3allSubj = nan(nsubj, nB,  ntypes, nLoc3);
    sep_3allSubj = nan(nsubj, nB,  nLoc3, ntypes);
    
    margORI_3allSubj = nan(nsubj, nB,  nLoc3, ntypes, nORI);
    margPredORI_3allSubj = nan(nsubj, nB,  nLoc3, ntypes, nORI);
    margR2ORI_3allSubj = nan(nsubj, nB,  nLoc3, ntypes);
    margParams_ORI_3allSubj = nan(nsubj, nB,  nLoc3, ntypes, length(paramInd_perF{1}));
    margTuningC_ORI_3allSubj = nan(nsubj, nB,  nLoc3, ntypes, nTuningC_ORI);
    
    margSF_3allSubj = nan(nsubj, nB,  nLoc3, ntypes, nSF);
    margPredSF_3allSubj = nan(nsubj, nB,  nLoc3, ntypes, nSF);
    margR2SF_3allSubj = nan(nsubj, nB,  nLoc3, ntypes);
    margParams_SF_3allSubj = nan(nsubj, nB,  nLoc3, ntypes, length(paramInd_perF{2}));
    margTuningC_SF_3allSubj = nan(nsubj, nB,  nLoc3, ntypes, nTuningC_SF);
    
    %%%%%%%%%
    for isubj = 1:nsubj
        subjName = subjList{isubj}; fprintf('  %s (%d/%d)...', subjName, isubj, nsubj)
        nblocks = nblocks_allSubj(isubj);
        
        for ibatch = 1:nBatches
            iB_start = (ibatch-1)*nB_perBatch+1;
            iB_end = ibatch*nB_perBatch;
            
            % Get data on HM
            load(sprintf('Data_OOD/%s%d/%s_batch%d_B%d_L67_N%d_%d_%d%s_ORI%d_SF%d.mat', ...
                subjName, nblocks, subjName, ibatch, nB_perBatch, flag_standEnergy, nORI, nSF, text_cut, ifamily_perF));
            iiLoc3 = 1; iiLoc2 = 1;
            %---------------%
            fxn_savePerLoc3
            %---------------%
            
            % Get data on LVM and UVM
            load(sprintf('Data_OOD/%s%d/%s_batch%d_B%d_L53_N%d_%d_%d%s_ORI%d_SF%d.mat', ...
                subjName, nblocks, subjName, ibatch, nB_perBatch, flag_standEnergy, nORI, nSF, text_cut, ifamily_perF));
            iiLoc3 = 2; iiLoc2 = 1;
            %---------------%
            fxn_savePerLoc3
            %---------------%
            
            iiLoc3 = 3; iiLoc2 = 2;
            %---------------%
            fxn_savePerLoc3
            %---------------%
            
        end
        fprintf('DONE\n')
    end % isubj
    % clear metrics_3allSubj
    
    fprintf('\nSaving...'), tic
    save(nameFileLoc3, '*_3allSubj')
else
    fprintf('\nLoading...'), tic
    load(nameFileLoc3, '*_3allSubj')
end
dur = toc; fprintf('  DONE (Dur %.1f min)\n', dur/60)
close all

%% caluclate coefficient of variation for critical measurements
% load(sprintf('%s/%s/n%d_n%d_Lall3_N_temp%d_B1_nA%d_conv1_IV1.mat', ...
%     nameFolderCompile_NOM, fileName_NOM, nsubj, ni, templateType, nModelsA))
% 
% getCoV % only when testing L6,5,3
% close all

%% similarity and separability
data_3allSubj = similarityAfterMirroring_3allSubj;
n = 'similarity';
% plot4_sep

data_3allSubj = sep_3allSubj;
n = 'separability';
plot4_sep

%% 1. Tuning function at three loc
close all

folderName  = sprintf('%s/%s/idvd/tuningFxn/Lall3', nameFigFolder, name_numFilters_Fitting); folderDir = dir(folderName); if isempty(folderDir), mkdir(folderName), end
sz_font = 10;

for ifeature=1:2
    if ifeature == 1
        ylim_tuning = [-.05, .2]; % left y-axis, limit
        ytext_R2 = [.19, .17, .15]; % left y-axis, the height of R2 for two loc
        marg_allSubj = margORI_3allSubj;
        margPred_allSubj = margPredORI_3allSubj;
        margR2_allSubj = margR2ORI_3allSubj;
    else
        ylim_tuning = [-.05, .15];
        ytext_R2 = [.14,  .12, .1]; % left y-axis, the height of R2 for two loc
        marg_allSubj = margSF_3allSubj;
        margPred_allSubj = margPredSF_3allSubj;
        margR2_allSubj = margR2SF_3allSubj;
    end
    
    xaxis = axis_tuning{ifeature};
    nfilters = length(xaxis);
    
    figure('Position', [0, 0, 2e3, 1.8e3])
    
    for isubj = 1:nsubj
        subjName = subjList{isubj};
        
        subplot(3,5, isubj), hold on
        
        %     marg_CI_2 = nan(nfilters, 2); % 1=lb of F/HM/LVM, 2=ub of P/VM/UVM
        
        for iiLoc = 1:nLoc3
            color = colors_comb(iLocComb_all(iiLoc), :);
            [marg_med, ~, ~, marg_CI_neg, marg_CI_pos] = getCI(marg_allSubj(isubj, :, iiLoc, itype, :), 1, 2, 1,1, CI_ratio);
            [margPred_med, margPred_lb, margPred_ub, ~, ~] = getCI(margPred_allSubj(isubj, :, iiLoc, itype, :), 1, 2, 1,1, CI_ratio);
            [~, margR2_lb, margR2_ub] = getCI(margR2_allSubj(isubj, :, iiLoc, itype), 1, 2, 1, 1, CI_ratio);
            
            %         if iiLoc == 1, marg_CI_2(:, iiLoc) = marg_CI_neg;
            %         else, marg_CI_2(:, iiLoc) = marg_CI_pos;
            %         end
            
            % extra lines
            yline(0, 'handlevisibility', 'off', 'linewidth', 1, 'color', [.7, .7, .7]);
            xline(ifeature-1, 'handlevisibility', 'off', 'linewidth', 1, 'color', [.7, .7, .7]);
            
            % plot kernels
            errorbar(xaxis, marg_med, marg_CI_neg, marg_CI_pos, 'color', color, 'CapSize',0, 'linestyle', 'none', 'linewidth', 2)
            %         plot(xaxis, marg_med, 'o', 'color', color,'MarkerEdgeColor', 'w' ,'linewidth', 1)
            
            % plot tuning fxn
            xaxis_itp = xaxis;
            plot(xaxis_itp, margPred_med', '-', 'color', color, 'linewidth', 2)
            %         patch([xaxis_itp, flip(xaxis_itp)], [margPred_lb', flip(margPred_ub')], color, 'FaceAlpha', .3, 'linestyle', 'none')
            
            % R2
            xtext_R2 = xaxis(18);
            text(xtext_R2, ytext_R2(iiLoc), sprintf('R^2 = [%d%%, %d%%]', round(margR2_lb*100), round(margR2_ub*100)), 'color', color)
            
        end % end of iiLoc
        
        % yaxis
        ylim(ylim_tuning)%, yticks(round(linspace(ylim_tuning(1), ylim_tuning(2), 5), 2))
        
        % xaxis
        xlim(axisLim{ifeature})
        xticks(axisTicks_tuning{ifeature})
        xticklabels(axisTL_tuning{ifeature})
        
        if isubj>1, xticklabels([]), yticklabels([]), end
        
        title(sprintf('S%d', isubj))
        
    end % isubj
    
    sgtitle(namesType{itype})
    set(findall(gcf, '-property', 'fontsize'), 'fontsize', sz_font)
    saveas(gcf, sprintf('%s/n%d_%s_%s.jpg', folderName, nsubj, namesFeature{ifeature}, namesType{itype}))
end

close all


%% 9a. [MAIN] corr between HVA/VMA: CS and tunC
clc, close all
iLocComb_all_all = {[6,7], [5,3]};
% iLocComb_all_all = {[6,7]};
for flag_noOL = 0%[0,1]
    for ii = iLocComb_all_all
        nameVarX = 'CS'; 
        nameVarY = 'tunC'; % 'pA', 'tunC', 'NOMparams'
        
        plot9_corrAsym
        
    end
end

%% 9b. [MAIN] corr between HVA/VMA: CS and pA
close all
iLocComb_all_all = {[6,7], [5,3]};
for flag_noOL = 0%[0,1]
    for ii = iLocComb_all_all
        nameVarX = 'CS'; % 'CS', 'pA'
        nameVarY = 'pA'; % 'pA', 'tunC', 'NOMparams'
        
        plot9_corrAsym
        
    end
end

%% [MAIN] 9c. corr between HVA/VMA: CS and NOMparams
close all
iLocComb_all_all = {[6,7], [5,3]};
for flag_noOL = 0%[0,1]
    for ii = iLocComb_all_all
        nameVarX = 'CS'; % 'CS', 'cpA'
        nameVarY = 'NOMparams'; % 'pA', 'tunC', 'NOMparams'
        
        plot9_corrAsym
    end
end

%% [MAIN] 9d. corr between HVA/VMA: ORI vs. SF
close all

folderNameFig_ORISF = sprintf('%s/%s/corrAsym/ORIvsSF', nameFigFolder, name_numFilters_Fitting); if isempty(dir(folderNameFig_ORISF)), mkdir(folderNameFig_ORISF), end

%%
iLocComb_all_all = {[6,7], [5,3]};

for ill=1:2
    iLocComb_all = iLocComb_all_all{ill};
    plot9d_ORIvsSF
end

%% 10. corr between HVA and VMA for CS/pA/tunC/NOMparams
close all
flag_noOL = 0;
nameVar = 'CS'; plot10_corr_HVA_VMA
nameVar = 'pA'; plot10_corr_HVA_VMA
nameVar = 'tunC_ORI'; plot10_corr_HVA_VMA
nameVar = 'tunC_SF'; plot10_corr_HVA_VMA
nameVar = 'NOMparams'; plot10_corr_HVA_VMA

%% 8a. [MAIN] correlation between NOM nLL and TunC
nameVarX = 'NOMnLL';
nameVarX = 'NOMnLLdelta';
nameVarY = 'tunC';

load('n13_NOM_nLL_L4', 'NOM_GoF_allSubj_allM_med')


for flag_zeroMean = 0
    plot8_corr
end

%% 8a. [MAIN] correlation between CS and TunC
nameVarX = 'CS';
nameVarY = 'tunC';

for flag_zeroMean = 0:2
    plot8_corr
end

%% 8b. [MAIN] correlation between CS and pA
nameVarX = 'CS';
nameVarY = 'pA';
for flag_zeroMean = 0:2
    plot8_corr
end

%% 8c. [MAIN] correlation between CS and estimated NOM params
nameVarX = 'CS';
nameVarY = 'NOMparams';
for flag_zeroMean = 0:2
    plot8_corr
end

%% 8d. [MAIN] correlation between pA and estimated NOM params
nameVarX = 'pA';
nameVarY = 'NOMparams';

for flag_zeroMean = 0:2
    plot8_corr
end

%% 8e. [SUPP] correlation between pA and tunC
nameVarX = 'pA';
nameVarY = 'tunC';

plot8_corr

%% 8f1. [SUPP] correlation between NOMparams (induced noise) and tunC
nameVarX = 'NOMparam1';
nameVarY = 'tunC';

plot8_corr

%% 8f2. [SUPP] correlation between NOMparams (constant noise) and tunC
nameVarX = 'NOMparam2';
nameVarY = 'tunC';

plot8_corr


%% LME
% to assess the contribution of sensitivity and internal noise to HVA/VMA
fprintf('\n\n*** ASYM of threshold predicted by that of orientation S (+), SF S (+), induced (-) and constant (-) internal noise [with random effect] ***\n')
getAsym = @(a) (a(:,1)-a(:,2))./(a(:,1)+a(:,2));
getAsym = @(a) (a(:,1)-a(:,2))./(a(:,2));
% getAsym = @(a) (a(:,1)-a(:,2))./(a(:,1));
iModelB = 2;
flag_asym = 1;
iLocComb = [6,7];% HVA
iLocComb = [5,3];% VMA

load(sprintf('%s/n%d_B1000_perf_L%d%d.mat', nameFolderCompile_BEHAV, nsubj, iLocComb))
CS_allSubj = getCI(cs_allSubj, 1, 2);
load(sprintf('%s/%s/n%d_B1000_L%d%d_%d_%d_m_ORI%dSF%d.mat', nameFolderCompile_RC, fileName, nsubj, iLocComb, nORI, nSF, ifamily_perF), 'margTuningC_ORI_allSubj', 'margTuningC_SF_allSubj')
ORI_S_allSubj = getCI(margTuningC_ORI_allSubj(:, :, :, itype, 2), 1, 2);
ORI_B_allSubj = getCI(margTuningC_ORI_allSubj(:, :, :, itype, 5), 1, 2);
SF_S_allSubj = getCI(margTuningC_SF_allSubj(:, :, :, itype, 2), 1, 2);
SF_B_allSubj = getCI(margTuningC_SF_allSubj(:, :, :, itype, 3), 1, 2);

load(sprintf('%s/%s/n%d_n%d_L%d%d_N_temp%d_B%d_nA%d_conv1_IV1.mat', ...
    nameFolderCompile_NOM, fileName_NOM, nsubj,  ni, iLocComb, templateType, iModelB, nModelsA))
switch iModelB
    case 1
        IndN_allSubj = getCI(squeeze(NOM_params_est_allSubj(:, :, iModelA, iModelB, :, 1)), 1, 3); %nsubj x nLoc x nModels x 1 x nB x nNOMparams
        ConN_allSubj = getCI(squeeze(NOM_params_est_allSubj(:, :, iModelA, iModelB, :, 2)), 1, 3); %nsubj x nLoc x nModels x 1 x nB x nNOMparams
    case 2
        ConN_allSubj = getCI(squeeze(NOM_params_est_allSubj(:, :, iModelA, iModelB, :, 1)), 1, 3); %nsubj x nLoc x nModels x 1 x nB x nNOMparam
end

if flag_asym
    tbl_asym = table(getAsym(CS_allSubj), getAsym(ORI_S_allSubj), getAsym(ORI_B_allSubj), getAsym(SF_S_allSubj), getAsym(ORI_B_allSubj), getAsym(IndN_allSubj), getAsym(ConN_allSubj), (1:nsubj)',...
        'VariableNames',{'CS_Asym','ORI_S_Asym', 'ORI_B_Asym', 'SF_S_Asym', 'SF_B_Asym', 'indN_Asym', 'conN_Asym', 'indSubj'});
    
    model_all = { {'ORI_S_Asym', 'SF_S_Asym'}, ... % just sensitivity
        {'ORI_S_Asym', 'ORI_B_Asym'}, ... % just ori
        {'SF_S_Asym', 'SF_B_Asym'}, ... % just SF
        {'ORI_S_Asym', 'SF_S_Asym', 'ORI_B_Asym', 'SF_B_Asym'}, ... % ori and SF, featural representation only
        {'indN_Asym', 'conN_Asym'}, ... % internal noise only
        {'ORI_S_Asym', 'SF_S_Asym', 'ORI_B_Asym', 'indN_Asym', 'conN_Asym'}, ... % all
        };
    
else
    loc = [ones(nsubj, 1), ones(nsubj, 1)*2];
    subj = [1:nsubj; 1:nsubj]';
    tbl = table(loc(:), CS_allSubj(:), ORI_S_allSubj(:), ORI_B_allSubj(:), SF_S_allSubj(:), ORI_B_allSubj(:),...
        IndN_allSubj(:), ConN_allSubj(:), subj(:),...
        'VariableNames',{'loc', 'CS','ORI_S', 'ORI_B', 'SF_S', 'SF_B', 'indN', 'conN', 'indSubj'});
    
    model_all = { {'ORI_S* loc', 'SF_S* loc'}, ... % just sensitivity
        {'ORI_S* loc', 'ORI_B* loc'}, ... % just ori
        {'SF_S* loc', 'SF_B* loc'}, ... % just SF
        {'ORI_S* loc', 'SF_S* loc', 'ORI_B* loc'}, ... % ori and SF, featural representation only
        {'indN* loc', 'conN* loc'}, ... % internal noise only
        {'ORI_S* loc', 'SF_S* loc', 'ORI_B* loc', 'indN* loc', 'conN* loc'}, ... % all
        };
end
nModels = length(model_all);

BIC = []; AIC=[];
for iModel = 1:nModels
    text_LME = [];
    for iVar = 1:length(model_all{iModel})
        text_LME = sprintf('%s %s +', text_LME, model_all{iModel}{iVar});
    end
    if flag_asym
        text_LME = sprintf('CS_Asym ~ %s (1|indSubj)', text_LME);
        lme = fitlme(tbl_asym, text_LME);
    else
        text_LME = sprintf('CS ~ %s (1|indSubj)', text_LME);
        lme = fitlme(tbl, text_LME);
    end
    
    BIC(iModel) = lme.ModelCriterion.BIC;
    AIC(iModel) = lme.ModelCriterion.AIC;
    
end
clc

BIC_delta = BIC - min(BIC)
AIC_delta = AIC- min(AIC)

iModel_opt = find(AIC_delta==0);
model_opt = model_all{iModel_opt};
nVar = length(model_opt);

text_LME = [];
for iVar = 1:nVar
    text_LME = sprintf('%s %s +', text_LME, model_opt{iVar});
end

if flag_asym
    text_LME = sprintf('CS_Asym ~ %s (1|indSubj)', text_LME);
    lme = fitlme(tbl_asym, text_LME);
else
    text_LME = sprintf('CS ~ %s (1|indSubj)', text_LME);
    lme = fitlme(tbl, text_LME);
end

text_symmary = [];
nItem = size(lme.Coefficients)-1;
for iItem = 1:nItem
    text_symmary = [text_symmary, sprintf('%s: %.2f [%.2f, %.2f], p=%.3f\n', ...
        lme.CoefficientNames{iItem+1}, lme.Coefficients.Estimate(iItem+1), lme.Coefficients.Lower(iItem+1), lme.Coefficients.Upper(iItem+1), lme.Coefficients.pValue(iItem+1))];
end

fprintf(text_symmary)

%%
close all
% end