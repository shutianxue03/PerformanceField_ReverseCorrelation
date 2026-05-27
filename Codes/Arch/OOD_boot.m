function OOD_boot(ibatch, isubj, iiLoc_all_all, nB_perBatch, ifamilyORI, ifamilySF, flag_cutMapping, flag_mirrorMapping)
clc
close all
warning off
format compact

addpath(genpath('fxn_analysis_RC_v2'))
addpath(genpath('fxn_MC'))
addpath(genpath('fxn_model'))
addpath(genpath('fxn_exp'))
addpath(genpath('Data_OOD'))

time_start = datetime('now')

%%
% nB=100: ~30 mins
% nB=1e3: 8-10h

%--------------%
SX_RC1_setting
%--------------%

subjList = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA', 'CS', 'DT', 'DU', 'RC', 'SR'};
nblocks_allSubj = [200, 240, 210, 240, 220, 220, 215, 220, 205, 210, 205, 205, 205, 195];

subjName = subjList{isubj};
nblocks = nblocks_allSubj(isubj);

iSess_start = 1; % sessions before this number are discardrd to obtain a high and stable quality of data
flag_PatchMode = 2; % 1=target patch; 2=noise patch
flag_standEnergy = 1; % 1=standardize energy; 0=do NOT

%% define the fitting function
ifamily_perF = [ifamilyORI, ifamilySF];
iLoc_all_all = [1,8; 6,7; 5,3];
iLocComb_all = iLoc_all_all(iiLoc_all_all, :);

paramInd_perF{1} = ones(1, length(namesParams_all{ifamilyORI})); % assuming all params free to vary
paramInd_perF{2} = ones(1, length(namesParams_all{ifamilySF})); % assuming all params free to vary

nTuningC_ORI = length(namesTunC_unit_perF{ifamilyORI, 2}); % preferred ori, peak amp, bandwidth, bottom
nTuningC_SF = length(namesTunC_unit_perF{ifamilySF, 2}); % peakSF, peak amp., bandwidth_full (octave), baseline

nLoc2 = length(iLocComb_all);
namesEnergySource = {'TARGET', 'NOISE'};

%% print
fprintf('\n=======================\n%s, nB=%d, nORI=%d x nSF=%d\n%s vs. %s\nEnergy derived from %s patch (standardization=%d)\nCut (1=cut edges): %d\nMirror (1-mirror mapping over 0 ORI): %d\n%s: Family#%d %s [%s]\n%s: Family#%d %s [%s]\n=======================\n\n', ...
    subjName, nB_perBatch, nORI, nSF, namesLocComb{iLocComb_all(1)}, namesLocComb{iLocComb_all(2)}, ...
    namesEnergySource{flag_PatchMode}, flag_standEnergy, ...
    flag_cutMapping, flag_mirrorMapping,...
    namesFeature{1}, ifamily_perF(1), namesFamily_all{ifamily_perF(1)}, num2str(paramInd_perF{1} ), ...
    namesFeature{2}, ifamily_perF(2), namesFamily_all{ifamily_perF(2)}, num2str(paramInd_perF{2} ))

%% get dir
nameFile_behav = sprintf('Data_OOD/%s%d/%s_behavMeas.mat', subjName, nblocks, subjName);
nameFile_energy = sprintf('Data_OOD/%s%d/%s_energy_N_%d_%d.mat', subjName, nblocks, subjName, nORI, nSF);

if flag_mirrorMapping, text_cut = '_m'; else, text_cut = ''; end
nameFile_Boot = sprintf('Data_OOD/%s%d/%s_batch%d_B%d_L%d%d_N%d_%d_%d%s_ORI%d_SF%d.mat', ...
    subjName, nblocks, subjName, ibatch, nB_perBatch, iLocComb_all, flag_standEnergy, nORI, nSF, text_cut, ifamily_perF);

%% load behavior data
fprintf('Loading behav data...'),tic
load(nameFile_behav, 'dataMatrix')
dur = toc; fprintf('DONE (Dur %.1f min)\n', dur/60)

%% load energy profiles
fprintf('Loading Source Energy...'), tic
load(nameFile_energy)
dur = toc; fprintf('DONE (Dur %.1f min)\n', dur/60)

if flag_PatchMode == 1, e3D_allT = e3D_target_allT; else, e3D_allT = e3D_noise_allT; end
[ntrials_allT, ~, ~] = size(e3D_allT);

nSess = ntrials_allT/500;
ntrialsAll = ntrials_allT/5; % ntrials per loc
ntrialsAll_s = (nSess-iSess_start+1)*100; 

%% empty holders
metrics_allB = nan(nB_perBatch, nLoc2, nmetrics); %(iB, iiLocComb, :) = [dprime, c, pC, pHit, pFA, pA3];
dataMtx_allB = nan(nB_perBatch, nLoc2, ntrialsAll_s*4, 11); % 4 because P has 4 single loc
pYES_tgt_allB = nan(nB_perBatch, nLoc2, ntypes, nbins_e);  %(iB, iiLocComb, :, :) = pYES_tgt;
ebin_tgt_allB = nan(nB_perBatch, nLoc2, ntypes, nbins_e); %(iB, iiLocComb, :, :) = ebin_tgt;
kernels2D_allB = nan(nB_perBatch, nLoc2, ntypes, nORI, nSF); %(iB, iiLocComb, :, :, :) = kernels2D;
sep_allB = nan(nB_perBatch, nLoc2, ntypes); %(iB, iiLocComb, :) = sep;
% secondary
pYES_tgt_norm_allB = nan(nB_perBatch, nLoc2, ntypes, nbins_e); %(iB, iiLocComb, :, :) = pYES_tgt_norm;
ebin_tgt_norm_allB = nan(nB_perBatch, nLoc2, ntypes, nbins_e); %(iB, iiLocComb, :, :) = ebin_tgt_norm;
intercept2D_allB = nan(nB_perBatch, nLoc2, ntypes, nORI, nSF); %(iB, iiLocComb, :, :, :) = intercept2D;
R2_2D_allB = nan(nB_perBatch, nLoc2, ntypes, nORI, nSF); %(iB, iiLocComb, :, :, :) = R2_2D;
R2_Tjur_allB = nan(nB_perBatch, nLoc2, ntypes, nORI, nSF); %(iB, iiLocComb, :, :, :) = R2_Tjur;
pValues_allB = nan(nB_perBatch, nLoc2, ntypes, nORI, nSF, 2); %(iB, iiLocComb, :, :, :, :) = pValues;
pCat_allB = nan(nB_perBatch, nLoc2, ntypes, nORI, nSF); %(iB, iiLocComb, :, :, :) = pCat;

similarityAfterMirroring_allB = nan(nB_perBatch, ntypes, nLoc2);
margORI_allB = nan(nB_perBatch, nLoc2, ntypes, nORI); %(iB, :, :, :) = margORI;
margPredORI_allB = nan(nB_perBatch, nLoc2, ntypes, nORI); %(iB, :, :, :) = margPred_ORI;
margParamsORI_allB = nan(nB_perBatch, nLoc2, ntypes, length(paramInd_perF{1})); %(iB, :, :, :) = margParams_ORI;
margTuningC_ORI_allB = nan(nB_perBatch,  nLoc2, ntypes, nTuningC_ORI);
margR2ORI_allB = nan(nB_perBatch, nLoc2, ntypes); %(iB, :, :, :) = margR2_ORI;

margSF_allB = nan(nB_perBatch, nLoc2, ntypes, nSF); %(iB, :, :, :) = margSF;
margPredSF_allB = nan(nB_perBatch, nLoc2, ntypes, nSF); %(iB, :, :, :) = margPred_SF;
margParamsSF_allB = nan(nB_perBatch, nLoc2, ntypes, length(paramInd_perF{2})); %(iB, :, :, :) = margParams_SF;
margTuningC_SF_allB = nan(nB_perBatch,  nLoc2, ntypes, nTuningC_SF);
margR2SF_allB = nan(nB_perBatch, nLoc2, ntypes); %(iB, :, :, :) = margR2_SF;

%%
for iB = 1:nB_perBatch
    fprintf('\n%d/%d...', iB, nB_perBatch) % no need to use 'parfor', which is used in fitting
    margORI_perComb = nan(nLoc2, ntypes, nORI);
    margSF_perComb = nan(nLoc2, ntypes, nSF);
    pref_ORI_perComb = nan(nLoc2, 1);
    
    for iiLocComb = 1:nLoc2
        iLocComb = iLocComb_all(iiLocComb);
        e3D_resampled = [];
        dataMtx_resampled = [];
        respC_resampled = [];
        
        % to equate the number of trials per loc, single location got doubled (HM vs. VM) even quandrupled (F vs. P)
        if iLocComb < 6, iLoc_all = ones(1,4)*iLocComb;
        else, switch iLocComb , case 6, iLoc_all = [2,2,4,4]; case 7, iLoc_all = [3,3,5,5]; case 8, iLoc_all = 2:5; end
        end
        
        nLoc_comb = length(iLoc_all);
        nChosen = ntrialsAll_s;
        iSess_select = iSess_start:nSess;
        
        for iiLoc = 1:nLoc_comb
            dataMtx_nonrand = [];
            e3D_nonrand = [];
            if iSess_select > 1 % selecting trials
                for iSS = 1:length(iSess_select)
                    indLoc = (dataMatrix(:, 5) == iLoc_all(iiLoc)) & (dataMatrix(:, 2) == iSess_select(iSS));
                    dataMtx_nonrand = [dataMtx_nonrand; dataMatrix(indLoc, :)];
                    e3D_nonrand = cat(1, e3D_nonrand, e3D_allT(indLoc, :, :));
                end % iSS
            else % without selecting trials
                indLoc = dataMatrix(:, 5) == iLoc_all(iiLoc);
                dataMtx_nonrand = [dataMtx_nonrand; dataMatrix(indLoc, :)];
                e3D_nonrand = cat(1, e3D_nonrand, e3D_allT(indLoc, :, :));
            end
            
            iPair_nonrand = dataMtx_nonrand(:, 8);
            iPair_unik_nonrand = unique(iPair_nonrand); % range: [1, ntrials_ofThisSubj/2]
            
            rng('shuffle')
            % get the itrial of selected data (for each loc) according to the unik pair
            indRand_unikPair = randi(ntrialsAll_s/2, [nChosen/2, 1]); % range: min>=1, max<=ntrialsAll/2
            [e3D_rand, dataMtx_rand, respC] = fxn_getRespC_simp(indRand_unikPair, iPair_unik_nonrand, e3D_nonrand, dataMtx_nonrand);
            
            e3D_resampled = [e3D_resampled; e3D_rand];
            dataMtx_resampled = [dataMtx_resampled; dataMtx_rand];
            respC_resampled = [respC_resampled; respC];
        end % iiLoc
        
        %% get behav matrix
        iPRS = dataMtx_resampled(:, 6);
        resp = dataMtx_resampled(:, 9);
        pHit = mean(iPRS & resp)*2;
        pFA = mean(~iPRS & resp)*2;
        pC = mean(iPRS == resp);
        [dprime, c] = SX_sim06_SDT(pHit, pFA);
        pA3 = nanmean(respC_resampled);
        
        %% standardize energy
        %------------------------------------------------------------------------------------------------%
        e3D_norm = normEnergy(e3D_resampled, dataMtx_resampled(:, 11), dataMtx_resampled(:, 6));
        %------------------------------------------------------------------------------------------------%
        if flag_standEnergy == 1, e3D = e3D_norm;
        else, e3D = e3D_resampled;
        end
        
        %% pYES vs. binned energy [SX_RC5_slope]
        % un-standardized energy
        [pYES_tgt, ebin_tgt] = SX_RC5_slope(dataMtx_resampled, e3D_resampled, ntypes, nbins_e);
        %         quickPlot5_slope
        
        % standardized energy
        [pYES_tgt_norm, ebin_tgt_norm] = SX_RC5_slope(dataMtx_resampled, e3D_norm, ntypes, nbins_e);
        %         fprintf(' Energy binned.')
        
        %% Probit regression to get KERNELs [SX_RC6_kernel]
        filtersSF_all = 2.^linspace(log2(noise.SF_low), log2(noise.SF_high), nSF);
        [kernels2D, intercept2D, margORI_, margSF_, R2_2D, R2_Tjur, pValues, pCat, sep] = ...% combine trials by averaging single loc
            SX_RC6_kernel_parfor(e3D, dataMtx_resampled, filtersSF_all, filtersOri_all);
        %         quickPlot6_kernel
        
        %% get preferred ORI (where kernels peak before mirroring)
        kernelORI = mean(kernels2D(2, :, :), 3);
        x_intp = linspace(axis_tuning{1}(1), axis_tuning{1}(end), 1e3);
        kernelORI_intp = interp1(axis_tuning{1}, kernelORI, x_intp);
        [~, imax] = max(kernelORI_intp);
        pref_ORI_perComb(iiLocComb) = x_intp(imax); 
        
        %% process kernels (cutting & mirroring)
        if flag_cutMapping, kernels2D = kernels2D(:, cut_ORI, cut_SF); end
        if flag_mirrorMapping
            indMir = (nORI-1)/2;
            kk_left = kernels2D(:, 1:indMir, :);
            kk_right = kernels2D(:, indMir+2:end, :);
            kk_mid = kernels2D(:, indMir+1,:);
            kk_ave = (kk_left + flip(kk_right,2))/2;
            kernels2D_mir = cat(2, kk_ave, kk_mid, flip(kk_ave, 2));
            
            % get similarity after mirroring
            for itype = 1:ntypes
                similarityAfterMirroring_allB(iB, itype, iiLocComb) = corr2(squeeze(kernels2D(itype, :, :)), squeeze(kernels2D_mir(itype, :, :)));
            end
            kernels2D = kernels2D_mir;
        end
        
        %% get margORI/SF_perComb  and calculate separability
        % NO need to save, as the reorganized version will be saved)
        for itype = 1:ntypes
            e2D_cut_ = squeeze(kernels2D(itype, :, : ));
            e_min = min(e2D_cut_(:));
            e2D_cut = e2D_cut_ - e_min + eps; % make all energy values positive
            margORI = mean(e2D_cut, 2);
            margSF = mean(e2D_cut, 1);
            e2D_recon = mtimes(margORI, margSF);
            sep(itype) = corr2(e2D_cut, e2D_recon);
            % marg
            margORI_perComb(iiLocComb, itype, :) = margORI + e_min; % stupid code, need to feed three types into SX_RC7_fitting
            margSF_perComb(iiLocComb, itype, :) = margSF + e_min;
        end % itype
        
        %% compile
        % primary
        metrics_allB(iB, iiLocComb, :) = [dprime, c, pC, pHit, pFA, pA3, nan]; % nan is for cst_ln
        dataMtx_allB(iB, iiLocComb, :, :) = dataMtx_resampled; % to extract RT and cst
        pYES_tgt_allB(iB, iiLocComb, :, :) = pYES_tgt;
        ebin_tgt_allB(iB, iiLocComb, :, :) = ebin_tgt;
        kernels2D_allB(iB, iiLocComb, :, :, :) = kernels2D;
        sep_allB(iB, iiLocComb, :) = sep;
        % secondary
        pYES_tgt_norm_allB(iB, iiLocComb, :, :) = pYES_tgt_norm;
        ebin_tgt_norm_allB(iB, iiLocComb, :, :) = ebin_tgt_norm;
        intercept2D_allB(iB, iiLocComb, :, :, :) = intercept2D;
        R2_2D_allB(iB, iiLocComb, :, :, :) = R2_2D;
        R2_Tjur_allB(iB, iiLocComb, :, :, :) = R2_Tjur;
        pValues_allB(iB, iiLocComb, :, :, :, :) = pValues;
        pCat_allB(iB, iiLocComb, :, :, :) = pCat;
        
    end % iiLocComb
    
    %% fit fxns to kernels [SX_RC7_fitting]
    flag_interpolate = 0;
    
    iF = 1; [margORI, margPred_ORI, margParams_ORI, margR2_ORI] = SX_RC7_fitting(iF, ...
        margORI_perComb, ub_full_all{ifamily_perF(iF)}, lb_full_all{ifamily_perF(iF)}, ...
        ifamily_perF, paramInd_perF, flag_standEnergy);
    iF = 2; [margSF, margPred_SF, margParams_SF, margR2_SF] = SX_RC7_fitting(iF, ...
        margSF_perComb, ub_full_all{ifamily_perF(iF)}, lb_full_all{ifamily_perF(iF)}, ...
        ifamily_perF, paramInd_perF, flag_standEnergy);
    
%     quickPlot7_fitting
    %     fprintf(' Kernel fitted.')
    
    %% compile
    margORI_allB(iB, :, :, :) = margORI;
    margPredORI_allB(iB, :, :, :) = margPred_ORI;
    margParamsORI_allB(iB, :, :, :) = margParams_ORI;
    margR2ORI_allB(iB, :, :, :) = margR2_ORI;
    
    margSF_allB(iB, :, :, :) = margSF;
    margPredSF_allB(iB, :, :, :) = margPred_SF;
    margParamsSF_allB(iB, :, :, :) = margParams_SF;
    margR2SF_allB(iB, :, :, :) = margR2_SF;
    
    %% extract tuning characteristics
    for iiLoc = 1:nLoc2
        for itype = 2%1:ntypes
            % ORI
            ifeature= 1;
            tuningC_ORI = fxn_getTuningC(axis_tuning{ifeature}, ifeature, ifamily_perF(ifeature), ...
                squeeze(margPred_ORI(iiLoc, itype,:)), squeeze(margParams_ORI(iiLoc, itype,:)));
            % SF
            ifeature=2;
            tuningC_SF = fxn_getTuningC(axis_tuning{ifeature}, ifeature, ifamily_perF(ifeature), ...
                squeeze(margPred_SF(iiLoc, itype,:)), squeeze(margParams_SF(iiLoc, itype,:)));
            % compile
            margTuningC_ORI_allB(iB, iiLoc, itype, :) = [pref_ORI_perComb(iiLoc), tuningC_ORI];
            margTuningC_SF_allB(iB, iiLoc, itype, :) = tuningC_SF;
            
        end % itype
    end % ii
    
end % end of iB

%% save
save(nameFile_Boot, '*_allB')

fprintf('ALL DONE\n')

%%
time_end = datetime('now')

time_end - time_start
