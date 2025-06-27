

%% preprocessing

% extract basic data and compute (1) behavioral measurements and (2) energy profile for each observer

clc
close all
warning off
format compact

addpath(genpath('Data'))
addpath(genpath('Data_OOD'))
addpath(genpath('fxn_exp'))
addpath(genpath('fxn_analysis_RC_v2'))
addpath(genpath('fxn_simulation'))

%% load params
%-------------%
SX_RC1_setting
%-------------%
subjList = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA',  'DT', 'CS', 'DU', 'SR', 'RC'};
nblocks_allSubj = [200, 240, 210, 240, 220, 220, 215, 220, 205, 210, 205, 205, 205, 195, 205];
nsubj = length(subjList);

ib_start = 1; % do NOT change!

%% create filters
[filter_sin, filter_cos] = SX_sim02_setFilters(stim, noise.filtersSF_all, filtersOri_all, 0);
fprintf('\nPool of filters created: nORI=%d, nSF=%d [%.2f, %.2f]\n', ...
    length(filtersOri_all), length(noise.filtersSF_all), noise.SF_low, noise.SF_high)

%%
for isubj = 1:nsubj
    subjName = subjList{isubj};
    nameFolder_Data = sprintf('Data/%s', subjName);
    nameFile_Data = sprintf('%s/%s_exp*', nameFolder_Data, subjName);
    dirFile_Data = dir(nameFile_Data); % find out all possible files
    nBlocks = length(dirFile_Data); % number of blocks finished
    
    %% if not consider all blocks
    ib_start = 1;
    if strcmp(subjName, 'SP'), ib_start = 81; end
    ib_end = nBlocks;
    nBlocks = ib_end - ib_start+1;
    fprintf('\n%s%d (S%d/%d) ...\n', subjName, nBlocks, isubj, nsubj)
    
    %% create file names and dir
    nameFolder_OOD =  sprintf('Data_OOD/%s%d/', subjName, nBlocks);
    if isempty(dir(nameFolder_OOD)), mkdir(nameFolder_OOD), end
    
    %% 1. extract basic data
    nameBehavMeas = sprintf('Data_OOD/%s%d/%s_behavMeas.mat', subjName, nBlocks, subjName);
    dirBehavMeas = dir(nameBehavMeas);
    nameTgtPatch = sprintf('Data/%s/%s_target.mat', subjName, subjName);
    nameNoisePatch = sprintf('Data/%s/%s_noise.mat', subjName, subjName);
    dirPatch_target = dir(nameTgtPatch);
    
    if isempty(dirBehavMeas)
        flagGetPatch = 1;
        %         if isempty(dirPatch_target), flagGetPatch = 1; else, flagGetPatch = 0;end
        % ======================
        SX_RC3a_extractSubjData
        SX_RC3b_organizeData_v2
        % ======================
    end
    
    % Since name* is deleted in SX_RC3b, define names below
    
    %% 2. calculate source energy (from Noise & target patch)
    %     the energy profile derived from noise patch, transformed or not, are NOT different
    name1 = {'N', 'T'};
    name2 = {'noise', 'target'};
    name3 = {'noise_allT', 'target_allT'};
    
    for iPatchMode = 2 % 1=noise patch; 2=target patch
        namePatch = sprintf('Data/%s/%s_%s.mat', subjName, subjName, name2{iPatchMode});
        nameEnergySource = sprintf('Data_OOD/%s%d/%s_energy_%s_%d_%d.mat', ...
            subjName, nBlocks, subjName, name1{iPatchMode}, nORI, nSF);
        
        % loading patches
        if ~exist(name3{iPatchMode}, 'var')
            tic
            fprintf('\n  Loading %s patches...', name2{iPatchMode})
            load(namePatch)
            dur = toc; fprintf('DONE (Dur %.1f min)\n', dur/60)
        end
        switch iPatchMode, case 1, patch_allT = noise_allT; case 2, patch_allT = target_allT; end
        assert(size(patch_allT, 1) == nBlocks*100)
        
        % generate energy profile
        fprintf('  Generating source energy (from %s) ...\n', name2{iPatchMode}), tic
        [e3D_allT, phase3D_allT] = SX_RC4_Energy_parfor(stim.mask, patch_allT, filter_sin, filter_cos);
        dur = toc; fprintf('  * Dur %.1f min *\n', dur/60)
        assert(size(e3D_allT, 1) == nBlocks*100)
        
        % saving
        fprintf('  Saving source energy (%s patch) ...', name2{iPatchMode}), tic
        switch iPatchMode
            case 1, e3D_noise_allT = e3D_allT; phase3D_noise_allT = phase3D_allT;
            case 2, e3D_target_allT = e3D_allT; phase3D_target_allT = phase3D_allT;
        end
        save(nameEnergySource, sprintf('e3D_%s', name3{iPatchMode}), sprintf('phase3D_%s', name3{iPatchMode}))
        dur = toc; fprintf('DONE (Dur %.1f min)\n', dur/60)
        clear *allT
    end % for iPatchMode
    
    %%
    clear patch* noise* % so that these two files won't be saved in the file of the next observer
    %     e3D = e3D_noise_allT; for iLoc = 1:5, quickPlot_energy, end
    %     e3D = e3D_target_allT; for iLoc = 1:5, quickPlot_energy, end
    
end % end of isubj

% plot contrast threshold, pA and accuracy as a function of session for each observer
plot_BehavMeasurebySess

