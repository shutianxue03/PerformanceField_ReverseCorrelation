% SX_analysis1_extractData.m
% This script extracts basic data from observers' data files and computes
% behavioral measurements and energy profile for each observer.
% It also sets up the necessary directories and parameters for further analysis.
% - It includes functions to extract data, organize it, and compute energy profiles.

% Analysis steps and the corresponding modular functions
% 1. SX_sim02_setFilters: Creates a bank of Gabor filters for deriving energy profiles.
% 2. SX_RC3a_extractSubjData: Extracts behavioral data and trial-wise patches from observers' data files.
% 3. SX_RC3b_organizeData_v2: Organizes the extracted data into a structured format.
% 4. SX_RC4_Energy_parfor: Computes the energy profile from noise and target patches.
% 5. plot_BehavMeasurebySess: Plots behavioral measurements and energy profiles.

% Outputs
% - Behavioral measurements and energy profiles saved in specified directories.
% - The script is designed to be run for each observer in the subject list.

%% Preprocessing
% extract basic data and compute (1) behavioral measurements and (2) energy profile for each observer

clear all, clc, close all, warning off, format compact

% addpath(genpath('Data'))
% addpath(genpath('Data_OOD'))
addpath(genpath('Codes/fxn_exp'))
addpath(genpath('Codes/fxn_analysis_RC_v2'))
addpath(genpath('Codes/fxn_simulation'))

%% load params
%-------------%
SX_RC1_setting
%-------------%
subjList = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA',  'DT', 'CS', 'DU', 'SR', 'RC'};
nblocks_allSubj = [200, 240, 210, 240, 220, 220, 215, 220, 205, 210, 205, 205, 205, 195, 205];
nsubj = length(subjList);

ib_start = 1; % do NOT change!

%% create filters
[filter_sin, filter_cos] = SX_sim02_setFilters(stim, noise.filtersSF_all, filtersOri_all, fxn_getSigma_SPdomain, 0);
fprintf('\nPool of filters created: nORI=%d, nSF=%d [%.2f, %.2f]\n', ...
    length(filtersOri_all), length(noise.filtersSF_all), noise.SF_low, noise.SF_high)

%%
for isubj = 1:nsubj
    subjName = subjList{isubj};
    nameFolder_DataData = sprintf('%s/Data/%s', nameFolder_Data, subjName);
    nameFile_Data = sprintf('%s/%s_exp*', nameFolder_DataData, subjName);
    dirFile_Data = dir(nameFile_Data); % find out all possible files
    nBlocks = length(dirFile_Data); % number of blocks finished
    
    %% if not consider all blocks
    ib_start = 1;
    if strcmp(subjName, 'SP'), ib_start = 81; end
    ib_end = nBlocks;
    nBlocks = ib_end - ib_start+1;
    fprintf('\n%s%d (S%d/%d) ...\n', subjName, nBlocks, isubj, nsubj)
    
    %% create file names and dir
    nameFolder_OOD =  sprintf('%s/Data_OOD_%d%d/%s%d/', nameFolder_Data, nORI, nSF, subjName, nBlocks);
    if isempty(dir(nameFolder_OOD)), mkdir(nameFolder_OOD), end
    
    %% 1. extract basic data
    nameBehavMeas = sprintf('%s/Data_OOD_%d%d/%s%d/%s_behavMeas.mat', ...
        nameFolder_Data, nORI, nSF, subjName, nBlocks, subjName);
    dirBehavMeas = dir(nameBehavMeas);
    nameFile_TgtPatch = sprintf('%s/Data/%s/%s_target.mat', nameFolder_Data, subjName, subjName);
    nameFile_NoisePatch = sprintf('%s/Data/%s/%s_noise.mat', nameFolder_Data, subjName, subjName);
    dirPatch_target = dir(nameFile_TgtPatch);
    
    if isempty(dirBehavMeas)
        % flag_getPatch = 1;
                if isempty(dirPatch_target), flag_getPatch = 1; else, flag_getPatch = 0;end
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
        nameFile_Patch = sprintf('%s/Data/%s/%s_%s.mat', nameFolder_Data, subjName, subjName, name2{iPatchMode});
        nameFile_EnergySource = sprintf('%s/Data_OOD_%d%d/%s%d/%s_energy_%s_%d_%d.mat', ...
            nameFolder_Data, nORI, nSF, subjName, nBlocks, subjName, name1{iPatchMode}, nORI, nSF);
        
        % loading patches
        if ~exist(name3{iPatchMode}, 'var')
            tic
            fprintf('\n  Loading %s patches...', name2{iPatchMode})
            load(nameFile_Patch)
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
        save(nameFile_EnergySource, sprintf('e3D_%s', name3{iPatchMode}), sprintf('phase3D_%s', name3{iPatchMode}))
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

