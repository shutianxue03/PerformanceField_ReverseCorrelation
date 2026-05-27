% Last updated on 07/14/2025 by Shutian Xue

% This script runs the RC analysis on all trials of each subject and derived kernels for paired locations: 
% Fovea vs. Periphery, HM vs. VM, Lower vs. Upper VM, and Left vs. Right HM
% The kernels are derived from the energy of the noise patches (or target patches)
% Figures are saved to the server: Figures/Kernels_2929

% Note: this script does not do resampling, so is OOD_boot

clc
close all
warning off
format compact

addpath(genpath('Codes/fxn_analysis_RC_v2'))
addpath(genpath('Codes/fxn_exp'))
addpath(genpath('Codes/fxn_MC'))

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
[filter_sin, filter_cos] = SX_sim02_setFilters(stim, noise.filtersSF_all, filtersOri_all, fxn_getSigma_SPdomain, 0);
fprintf('\n\n   Pool of filters created: nORI=%d, nSF=%d\n\n\n', length(filtersOri_all), length(noise.filtersSF_all))

%% MAIN LOOP
for isubj = 1:nsubj % Each subject takes <2min to run
    % get subject name and number of blocks
    subjName = subjList{isubj};
    nblocks = nblocks_allSubj(isubj);

    fprintf('%s%d (%d/%d) ...\n', subjName, nblocks, isubj, nsubj)

    % Define file names
    nameFile_BehavMeas = sprintf('%s/Data_OOD_%d%d/%s%d/%s_behavMeas.mat', nameFolder_Data, nORI, nSF, subjName, nblocks, subjName);
    nameDir_BehavMeas = dir(nameFile_BehavMeas);
    if flag_PatchMode == 1, namePatchMode = 'T'; else, namePatchMode = 'N'; end
    nameFile_SourceEnergy = sprintf('%s/Data_OOD_%d%d/%s%d/%s_energy_%s_%d_%d.mat', nameFolder_Data, nORI, nSF, subjName, nblocks, subjName, namePatchMode, nORI, nSF);
    nameFile_Kernel = sprintf('%s/Data_OOD_%d%d/%s%d/%s_kernel_%s_%d_%d.mat', nameFolder_Data, nORI, nSF, subjName, nblocks, subjName, namePatchMode, nORI, nSF);

    % Extract behav data matrix
    fprintf('  Loading Data matrix...'), tic
    load(nameFile_BehavMeas, 'dataMatrix')
    dur = toc; fprintf('DONE (Dur %.1f min)\n', dur/60)

    % Load the source energy (from which energy will be sampled in each bootstrapping)
    fprintf('  Loading Source Energy...'), tic
    load(nameFile_SourceEnergy)
    if flag_PatchMode == 1, e3D_allT = e3D_target_allT; else, e3D_allT = e3D_noise_allT; end
    dur = toc; fprintf('DONE (Dur %.1f min)\n', dur/60)

    % Check the data matrix
    [nTrials_allT, nORI, nSF] = size(e3D_allT);
    nTrialsPerLoc_original = nTrials_allT/5; % ntrials per loc
    nSess = nTrials_allT/500;

    % quick plot
    % e3D = e3D_target_allT; for iLoc = 1:5, quickPlot_energy, end

    fprintf('  Generating kernels per loc...   '), tic

    % Prepare variables
    pYES_tgt_perComb = cell(nLoc8, 1);
    ebin_tgt_perComb = pYES_tgt_perComb;
    pYES_tgt_norm_perComb = pYES_tgt_perComb;
    ebin_tgt_norm_perComb = pYES_tgt_perComb;
    kernels2D_perComb = pYES_tgt_perComb;
    intercept2D_perComb = pYES_tgt_perComb;
    margORI_perComb = pYES_tgt_perComb;
    margSF_perComb = pYES_tgt_perComb;

    % Loop through each loc combination
    for iLocComb = iLocComb_all
        if iLocComb < 6, iLoc_all = iLocComb;
        else, switch iLocComb , case 6, iLoc_all = [2,4]; case 7, iLoc_all = [3,5]; case 8, iLoc_all = 2:5; end
        end
        nLoc2 = length(iLoc_all);

        % Preallocate variables
        dataMtx_resampled = [];
        e3D_resampled = [];
        iPRS = [];
        cst = [];

        % Loop through each loc within the combination
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

        % Standardize energy
        %----------------------------------------------------%
        %         e3D_norm = normEnergy(e3D_resampled, iPRS, cst); % old
        e3D_norm = normEnergy(e3D_resampled, cst, iPRS); % corrected
        %----------------------------------------------------%

        % Derive pYES vs. binned energy [SX_RC5_slope]
        % With the un-standardized energy
        [pYES_tgt, ebin_tgt] = SX_RC5_slope(dataMtx_resampled, e3D_resampled, nTypes, nbins_e);

        % With the standardized energy
        [pYES_tgt_norm, ebin_tgt_norm] = SX_RC5_slope(dataMtx_resampled, e3D_norm, nTypes, nbins_e);

        % Conduct probit regression to derive KERNELs [SX_RC6_kernel]
        [kernels2D, intercept2D, margORI, margSF, R2_2D, R2_Tjur, pValues, pCat, sep] = ...% combine trials by averaging single loc
            SX_RC6_kernel_parfor(e3D_norm, dataMtx_resampled, filtersSF_all, filtersOri_all);

        % Compile results for each loc combination
        pYES_tgt_perComb{iLocComb} = pYES_tgt;
        ebin_tgt_perComb{iLocComb} = ebin_tgt;
        pYES_tgt_norm_perComb{iLocComb} = pYES_tgt_norm;
        ebin_tgt_norm_perComb{iLocComb} = ebin_tgt_norm;
        kernels2D_perComb{iLocComb} = kernels2D;
        intercept2D_perComb{iLocComb} = intercept2D;
        margORI_perComb{iLocComb} = margORI;
        margSF_perComb{iLocComb} = margSF;

        fprintf(' | %s ', namesLocComb{iLocComb})
    end % end of iLocComb

    dur = toc; fprintf('DONE (Dur %.1f min)\n', dur/60)

    % Save results
    fprintf('  Saving kernels...'), tic
    % save(nameFile_Kernel, '*_perComb', 'ntrialsPerLoc_*')
    save(nameFile_Kernel, '*_perComb') % deleted ntrialsPerLoc_* because it does not exist
    dur = toc; fprintf('DONE (Dur %.1f min)\n', dur/60)

end % end of isubj

%% Plot kernels of three loc pairs for each subject
close all, clc
namePatchMode = 'N'; % N=energy are computed from noise patches; T=energy are computed from target patches
iLocPair_all = [1, 8; 6, 7; 5, 3; 2, 4]; % F=1 vs. P=8; HM=6 vs. VM=7; LVM=5 vs. UVM=3; LHM=2 vs. RHM=4
[nLocPairs, nLoc2] = size(iLocPair_all);

% Loop through mirroring kernels or not
for flag_mirror = 0:1
    % Define whether to mirror the kernels
    if flag_mirror == 0, str_mirror = ''; else, str_mirror = '_mirrored'; end

    % Loop through each subject
    for isubj = 1:nsubj
        subjName = subjList{isubj};
        nblocks = nblocks_allSubj(isubj);

        % Load kernels
        nameFile_Kernel = sprintf('%s/Data_OOD_%d%d/%s%d/%s_kernel_%s_%d_%d.mat', nameFolder_Data, nORI, nSF, subjName, nblocks, subjName, namePatchMode, nORI, nSF);
        load(nameFile_Kernel, 'kernels2D_perComb')

        % Loop through each feature (1=ORI, 2=SF)
        for iFeature = 1:2

            % Define name of the figure folder
            nameFolder_Fig_kernel = sprintf('%s/Kernels_%d%d/%s', nameFolder_Figures, nORI, nSF, namesFeature{iFeature});
            if isempty(dir(nameFolder_Fig_kernel)), mkdir(nameFolder_Fig_kernel), end

            if iFeature == 1 % ORI
                figure('Position', [0 200 nLocPairs*250 600])
                marg = margORI_perComb;
                y_lim = [-.1, .25];
            else % SF
                figure('Position', [800 200 nLocPairs*250 600])
                marg  = margSF_perComb;
                y_lim = [-.05, .1];
            end

            isubplot = 1;
            for iType = 1:nTypes % 1=PRS, 2=ABS, 3=ALL trials
                % Loop through each loc pair
                for iiLocPair = 1:nLocPairs
                    iLocPair = iLocPair_all(iiLocPair, :);

                    % Mirror kernels if necessary
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

                    % Extract marginal tuning functions
                    marg1 = marg(iType, :, 1);
                    marg2 = marg(iType, :, 2);

                    % Plot the tuning functions of two loc
                    subplot(nTypes, nLocPairs, isubplot), hold on
                    plot(axis_tuning{iFeature},  marg1, '-', 'color', colors_comb(iLocPair(1), :), 'linewidth', 2)
                    plot(axis_tuning{iFeature},  marg2, '-', 'color', colors_comb(iLocPair(2), :), 'linewidth', 2)
                    ylim(y_lim), xline(iFeature-1, 'k-'); yline(0, 'k-');

                    xticks(axisTicks_tuning{iFeature})
                    xticklabels(axisTL_tuning{iFeature})
                    title(sprintf('[%s] %s vs. %s', namesType{iType}, namesLocComb{iLocPair(1)}, namesLocComb{iLocPair(2)}))

                    isubplot = isubplot+1;
                end % end of iiLocPair
            end % end of iType
            %         sgtitle(sprintf('[%s-%s] %s - %d/%d trials', namesFeature{iFeature}, namePatchMode, subjName, ntrialsPerLoc_selected, nTrialsPerLoc_original))
            sgtitle(sprintf('[%s-%s] %s (%s)', namesFeature{iFeature}, namePatchMode, subjName, str_mirror))

            % Save the figure (to the server)
            saveas(gcf, sprintf('%s/kernels_%s_%s%s.jpg', nameFolder_Fig_kernel,  subjName, namePatchMode, str_mirror))
        end % iFeature
    end % isubj

    close all
end % end of flag_mirror

fprintf('\n\n\n   All done! Kernels saved to %s\n\n', nameFolder_Fig_kernel)