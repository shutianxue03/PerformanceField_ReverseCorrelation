
% Step 2: [Local] run PF_RC_step2 to do RC analysis on ALL trials (no resampling)
% code is similar to OOD_boot (noresampling)
% use the sensitivity kernels derived from ALL trials to determine which model best captures ORI and SF

clc
close all
warning off
format compact

addpath(genpath('fxn_analysis_RC_v2'))
addpath(genpath('fxn_exp'))
addpath(genpath('fxn_MC'))
addpath(genpath('Data_OOD'))

%%
%---------------%
SX_RC1_setting
%---------------%
subjList = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA',  'DT', 'CS', 'DU', 'SR', 'RC'};
nblocks_allSubj = [200, 240, 210, 240, 220, 220, 215, 220, 205, 210, 205, 205, 205, 195, 205];
nSess = nblocks_allSubj/5;
nsubj = length(subjList);

iSess_start = 1; iSess_select = iSess_start:nSess;
flag_PatchMode = 2; % 1=use energy of targets; 2=energy of noise parts
iLocComb_all = 1:8;
fprintf('Discard data before # %d\n', iSess_start)

%% create gabor filters
[filter_sin, filter_cos] = SX_sim02_setFilters(stim, noise.filtersSF_all, filtersOri_all, 0);
fprintf('\n\n   Pool of filters created: nORI=%d, nSF=%d\n\n\n', length(filtersOri_all), length(noise.filtersSF_all))

%% MAIN LOOP
for isubj = 1:nsubj
    
    subjName = subjList{isubj};
    nblocks = nblocks_allSubj(isubj);
    
    fprintf('%s%d (%d/%d) ...\n', subjName, nblocks, isubj, nsubj)
    
    % get dir
    nameBehavMeas = sprintf('Data_OOD/%s%d/%s_behavMeas.mat', subjName, nblocks, subjName);
    dirBehavMeas = dir(nameBehavMeas);
    if flag_PatchMode == 1, namePatchMode = 'T'; else, namePatchMode = 'N'; end
    nameSourceEnergy = sprintf('Data_OOD/%s%d/%s_energy_%s_%d_%d.mat', subjName, nblocks, subjName, namePatchMode, nORI, nSF);
    nameKernel = sprintf('Data_OOD/%s%d/%s_kernel_%s_%d_%d.mat', subjName, nblocks, subjName, namePatchMode, nORI, nSF);
    
    % extract basic data matrix
    fprintf('  Loading Data matrix...'), tic
    load(nameBehavMeas, 'dataMatrix')
    dur = toc; fprintf('DONE (Dur %.1f min)\n', dur/60)
    
    % load the source energy (from which energy will be sampled in each bootstrapping)
    fprintf('  Loading Source Energy...'), tic
    load(nameSourceEnergy)
    if flag_PatchMode == 1, e3D_allT = e3D_target_allT; else, e3D_allT = e3D_noise_allT; end
    dur = toc; fprintf('DONE (Dur %.1f min)\n', dur/60)
    
    [nTrials_allT, nORI, nSF] = size(e3D_allT);
    nTrialsPerLoc_original = nTrials_allT/5; % ntrials per loc
    nSess = nTrials_allT/500;
    
    % quick plot
    % e3D = e3D_target_allT; for iLoc = 1:5, quickPlot_energy, end
    
    fprintf('  Generating kernels per loc...   '), tic
    % empty containers
    pYES_tgt_perComb = cell(nLoc8, 1);
    ebin_tgt_perComb = pYES_tgt_perComb;
    pYES_tgt_norm_perComb = pYES_tgt_perComb;
    ebin_tgt_norm_perComb = pYES_tgt_perComb;
    kernels2D_perComb = pYES_tgt_perComb;
    intercept2D_perComb = pYES_tgt_perComb;
    margORI_perComb = pYES_tgt_perComb;
    margSF_perComb = pYES_tgt_perComb;
    
    for iLocComb = iLocComb_all
        if iLocComb < 6, iLoc_all = iLocComb;
        else, switch iLocComb , case 6, iLoc_all = [2,4]; case 7, iLoc_all = [3,5]; case 8, iLoc_all = 2:5; end
        end
        nLoc2 = length(iLoc_all);
        
        dataMtx_resampled = [];
        e3D_resampled = [];
        iPRS = [];
        cst = [];
        
        for iiLoc = 1:nLoc2
            if iSess_start > 2 % selecting trials
                for iSS = 1:length(iSess_select)
                    indLoc = (dataMatrix(:, 5) == iLoc_all(iiLoc)) & (dataMatrix(:, 2) == iSess_select(iSS));
                    dataMtx_resampled = [dataMtx_resampled; dataMatrix(indLoc, :)];
                    e3D_resampled = cat(1, e3D_resampled, e3D_allT(indLoc, :, :));
                    iPRS = [iPRS; dataMatrix(indLoc, 6)];
                    cst = [cst; dataMatrix(indLoc, 11)];
                end % iSS
            else % without selecting trials
                indLoc = dataMatrix(:, 5) == iLoc_all(iiLoc);
                dataMtx_resampled = [dataMtx_resampled; dataMatrix(indLoc, :)];
                e3D_resampled = cat(1, e3D_resampled, e3D_allT(indLoc, :, :));
                iPRS = [iPRS; dataMatrix(indLoc, 6)];
                cst = [cst; dataMatrix(indLoc, 11)];
            end
            assert(length(unique([size(dataMtx_resampled, 1), size(e3D_resampled, 1), size(iPRS, 1), size(cst , 1)]))==1)
        end % iiLoc
        
        % standardize energy
        %----------------------------------------------------%
        %         e3D_norm = normEnergy(e3D_resampled, iPRS, cst); % old
        e3D_norm = normEnergy(e3D_resampled, cst, iPRS); % corrected
        %----------------------------------------------------%
        
        % pYES vs. binned energy [SX_RC5_slope]
        % un-standardized energy
        [pYES_tgt, ebin_tgt] =  SX_RC5_slope(dataMtx_resampled, e3D_resampled, nTypes, nbins_e);
        
        % standardized energy
        [pYES_tgt_norm, ebin_tgt_norm] = SX_RC5_slope(dataMtx_resampled, e3D_norm, nTypes, nbins_e);
        
        % Probit regression to get KERNELs [SX_RC6_kernel]
        [kernels2D, intercept2D, margORI, margSF, R2_2D, R2_Tjur, pValues, pCat, sep] = ...% combine trials by averaging single loc
            SX_RC6_kernel_parfor(e3D_norm, dataMtx_resampled, filtersSF_all, filtersOri_all);
        
        % compile vars for each iB
        pYES_tgt_perComb{iLocComb} = pYES_tgt;
        ebin_tgt_perComb{iLocComb} = ebin_tgt;
        pYES_tgt_norm_perComb{iLocComb} = pYES_tgt_norm;
        ebin_tgt_norm_perComb{iLocComb} = ebin_tgt_norm;
        kernels2D_perComb{iLocComb} = kernels2D;
        intercept2D_perComb{iLocComb} = intercept2D;
        margORI_perComb{iLocComb} = margORI;
        margSF_perComb{iLocComb} = margSF;
        
        fprintf('%s ', namesLocComb{iLocComb})
    end % iLocComb
    
    dur = toc; fprintf('DONE (Dur %.1f min)\n', dur/60)
    
    % save kernels to OOD
    fprintf('  Saving kernels...'), tic
    save(nameKernel, '*_perComb', 'ntrialsPerLoc_*')
    dur = toc; fprintf('DONE (Dur %.1f min)\n', dur/60)
    
end % isubj

%% plot kernels of three loc pairs
close all, clc
flag_mirror = 0;
namePatchMode = 'N';
iLocPair_all = [1, 8; 6, 7; 5, 3; 2, 4]; % F vs. P; HM vs. VM; LVM vs. UVM; LHM vs. RHM
[nLocPairs, nLoc2] = size(iLocPair_all);

for isubj = 1:nsubj
    subjName = subjList{isubj};
    nblocks = nblocks_allSubj(isubj);
    nameKernel = sprintf('Data_OOD/%s%d/%s_kernel_%s_%d_%d.mat', subjName, nblocks, subjName, namePatchMode, nORI, nSF);
    load(nameKernel, 'kernels2D_perComb')
    
    for iFeature = 1:2
        if iFeature == 1
            figure('Position', [0 200 nLocPairs*250 600])
            marg = margORI_perComb;
            y_lim = [-.1, .25];
        else
            figure('Position', [800 200 nLocPairs*250 600])
            marg  = margSF_perComb;
            y_lim = [-.05, .1];
        end
        
        isubplot = 1;
        for iType = 1:nTypes
            for iiLocPair = 1:nLocPairs
                iLocPair = iLocPair_all(iiLocPair, :);
                
                % mirror
                kernels2D_plot = cat(4, kernels2D_perComb{iLocPair(1)}, kernels2D_perComb{iLocPair(2)}); % nTypes x nORI x nSF x 2
                if flag_mirror==1
                    kernels2D_mir = nan(size(kernels2D_plot));
                    indMir = (nORI-1)/2;
                    kk_left = kernels2D_plot(:, 1:indMir, :, :);
                    kk_right = kernels2D_plot(:, indMir+2:end, :, :);
                    kk_mid = kernels2D_plot(:, indMir+1,:, :);
                    kk_ave = (kk_left + flip(kk_right, 2))/2;
                    kernels2D_plot = cat(2, kk_ave, kk_mid, flip(kk_ave,2));
                end
                switch iFeature
                    case 1, marg = squeeze(mean(kernels2D_plot, 3));
                    case 2, marg = squeeze(mean(kernels2D_plot, 2));
                end
                
                marg1 = marg(iType, :, 1);
                marg2 = marg(iType, :, 2);
                
                % tuning functions of two loc
                subplot(nTypes, nLocPairs, isubplot), hold on
                plot(axis_tuning{iFeature},  marg1, '-', 'color', colors_comb(iLocPair(1), :), 'linewidth', 2)
                plot(axis_tuning{iFeature},  marg2, '-', 'color', colors_comb(iLocPair(2), :), 'linewidth', 2)
                ylim(y_lim), xline(iFeature-1, 'k-'); yline(0, 'k-');
                
                xticks(axisTicks_tuning{iFeature})
                xticklabels(axisTL_tuning{iFeature})
                title(sprintf('[%s] %s vs. %s', namesType{iType}, namesLocComb{iLocPair(1)}, namesLocComb{iLocPair(2)}))
                
                isubplot = isubplot+1;
            end % iiLocPair
        end % iType
        %         sgtitle(sprintf('[%s-%s] %s - %d/%d trials', namesFeature{iFeature}, namePatchMode, subjName, ntrialsPerLoc_selected, nTrialsPerLoc_original))
        sgtitle(sprintf('[%s-%s] %s', namesFeature{iFeature}, namePatchMode, subjName))
        nameFolder_kernel = sprintf('Fig_kernels/%s', namesFeature{iFeature});
        if isempty(dir(nameFolder_kernel)), mkdir(nameFolder_kernel), end
        saveas(gcf, sprintf('%s/kernels_%s_%s.jpg', nameFolder_kernel,  subjName, namePatchMode))
    end % iFeature
end % isubj

close all