
% SX_analysis4_Boot.m
% Last updated by Shutian Xue on 07/21/2025
%
% Description:
%   This script compiles and analyzes bootstrapped data for each subject.
%   For each observer, it loads bootstrapped data (run on OOD), compiles individual and group-level metrics,
%   and saves the results for downstream analysis and plotting.
%
% Usage:
%   Run this script after bootstrapping (OOD_boot_current) is completed

clc, close all, clear all, warning off, format compact

addpath(genpath('fxn_MC'))
addpath(genpath('fxn_exp'))
addpath(genpath('fxn_analysis_RC_v2'))
addpath(genpath('fxn_simulation'))
addpath(genpath('fxn_RCplot'))
addpath(genpath('XueLandyCarrasco'))

%% basic settings
subjList = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA',  'DT', 'CS', 'DU', 'SR', 'RC'};
nblocks_allSubj = [200, 240, 210, 240, 220, 220, 215, 220, 205, 210, 205, 205, 205, 195, 205];

nsubj = length(subjList);

global flag_cutMapping ifamily_perF

nModelsA = 5; %CHECK THIS; SHOULD BE 3

%---------------%
SX_RC1_setting
%---------------%
flag_ABSprefORI = 1;%input('     >>> take abs of pref ORI? (1=YES, 0=NO): '); % 1=take abs, meanig distance from the signal; 0=the sign (+ or -) indicates bias
flag_cutMapping = 0; % 1=load original data (ORI29 SF29) and cut the edge
flag_mirrorMapping = 1; % 1=mirror 2D mapping over 0 deg ori
nB_perBatch = 1e2; % MUST be the same as OOD_boot_current.m!!
nBatches = 10; % MUST be the same as OOD_boot_current.m!!
nB = nB_perBatch * nBatches;
ni = 1e3;
flag_PatchMode = 2; % energy derives from target; 2=noise patch % MUST be the same as shell_all_OOD.sh!!
flag_block200 = 0; % 1; all observers are forced to have 200 blocks; 0=no
flag_standEnergy = 1; if flag_PatchMode==1, flag_standEnergy =1; end

% Determine the list of subjects
nsubj = length(subjList);
namesEnergySource = {'TARGET', 'NOISE'};
if flag_PatchMode == 1, namePatchMode = 'T';
else, namePatchMode = sprintf('N%d', flag_standEnergy);
end

if flag_mirrorMapping, text_cut = '_m'; else, text_cut = ''; end

%% define tuning fxn fitting params
% edit plotAll_part2.m % saved in XueLandyCarrasco
% edit plotAll_part1.m % saved in XueLandyCarrasco
%---------------%
SX_RC1_setting
%---------------%
fprintf('\n\nnORI=%d nSF=%d\n\n', nORI, nSF)

iLocComb_all_all = {[6,7], [5,3]};
% [1,8]: fovea vs. periphery
% [2,4]: LHM vs. RHM (also plot D vs. nonD) % run this before [6,7]!!
% [6,7]: HM vs. VM (also plot HM (D) vs. VM)
% [5,3]: LVM vs. UVM

ill =2;% input('     >>> Which location to compare (1=L18, 2=L67, 3=L53, 4=L24, 5=L47): ');
iLocComb_all = iLocComb_all_all{ill};% delete

ifamily_perF = [1,2]; %input('     >>> What is the model family for ORI and SF (e.g., [1,2]): '); % ORI: scaled Gaussian; SF: double peak

nTuningC_ORI = length(namesTunC_unit_perF{ifamily_perF(1), 2}); % preferred ori, peak amp, bandwidth, bottom, // pred_half
nTuningC_SF = length(namesTunC_unit_perF{ifamily_perF(2), 2}); % peakSF, peak amp., bandwidth_full (octave), baseline // bandwidth_L, bandwidth_R, bottom, pred_half, full bandwidth (linear scale)

paramInd_perF{1} = ones(1,length(namesParams_all{ifamily_perF(1)})); %
paramInd_perF{2} = ones(1,length(namesParams_all{ifamily_perF(2)}));

% Print summary of settings
fprintf('\nnsubj = %d, nB = %d\n\n%s vs. %s\nEnergy derived from %s patch (standardization=%d)\n\nORI: M%d %s [%s]\nSF: M%d %s [%s]\n', ...
    nsubj, nB, ...
    namesLocComb{iLocComb_all(1)}, namesLocComb{iLocComb_all(2)}, ...
    namesEnergySource{flag_PatchMode}, flag_standEnergy, ...
    ifamily_perF(1), namesFamily_all{ifamily_perF(1)}, num2str(paramInd_perF{1}), ...
    ifamily_perF(2), namesFamily_all{ifamily_perF(2)}, num2str(paramInd_perF{2}))

%% Define directories
nameFolder_compile = sprintf('%s/Data_compile/RC_%d%d%s', nameFolder_Data, nORI, nSF, text_cut);
nameFile_bootGroup = sprintf('%s/n%d_B%d_L%d%d_%d_%d%s_ORI%dSF%d.mat', ...
    nameFolder_compile, nsubj, nB, iLocComb_all, nORI, nSF, text_cut, ifamily_perF(1), ifamily_perF(2));

nameDir_bootGroup = dir(nameFile_bootGroup);

if isempty(nameDir_bootGroup)

    % Create directory
    mkdir(nameFolder_compile);

    % Preallocate variables
    %     if nmetrics == 9, nmetrics=8; end
    %     metrics_allSubj = nan(nsubj, nB,  nLoc2, nmetrics);
    %     cst_allSubj = nan(nsubj, nB,  nLoc2);
    %     RT_allSubj = nan(nsubj, nB,  nLoc2);
    pYES_tgt_allSubj = nan(nsubj, nB,  nLoc2, nTypes, nbins_e);
    ebin_tgt_allSubj = nan(nsubj, nB,  nLoc2, nTypes, nbins_e);
    similarityAfterMirroring_allSubj = nan(nsubj, nB,  nTypes, nLoc2);
    kernels2D_allSubj = nan(nsubj, nB,  nLoc2, nTypes, nORI, nSF);
    sep_allSubj = nan(nsubj, nB,  nLoc2, nTypes);

    margORI_allSubj = nan(nsubj, nB,  nLoc2, nTypes, nORI);
    margORIraw_allSubj = nan(nsubj, nB,  nLoc2, nTypes, nORI);
    margPredORI_allSubj = nan(nsubj, nB,  nLoc2, nTypes, nORI);
    margParams_ORI_allSubj = nan(nsubj, nB,  nLoc2, nTypes, length(paramInd_perF{1}));
    margTuningC_ORI_allSubj = nan(nsubj, nB,  nLoc2, nTypes, nTuningC_ORI);
    margR2ORI_allSubj = nan(nsubj, nB,  nLoc2, nTypes);

    margSF_allSubj = nan(nsubj, nB,  nLoc2, nTypes, nSF);
    margPredSF_allSubj = nan(nsubj, nB,  nLoc2, nTypes, nSF);
    margParams_SF_allSubj = nan(nsubj, nB,  nLoc2, nTypes, length(paramInd_perF{2}));
    margTuningC_SF_allSubj = nan(nsubj, nB,  nLoc2, nTypes, nTuningC_SF);
    margR2SF_allSubj = nan(nsubj, nB,  nLoc2, nTypes);

    % Loop through subjects
    for isubj = 1:nsubj

        % Extract and print subject info
        subjName = subjList{isubj};
        nblocks = nblocks_allSubj(isubj);
        fprintf('\n%s (%d/%d)...\n', subjName, isubj, nsubj)

        % Load behav metric of each subject
        nameFile_behavMeas = sprintf('%s/%s%d/%s_behavMeas.mat', nameFolder_Data_OOD, subjName, nblocks, subjName);
        fprintf('Loading behav data...')
        load(nameFile_behavMeas)
        fprintf(' DONE\n')

        % Load bootstrapped outputs of each subject
        if flag_cutMapping, nORI_load = nan; nSF_load = nan;
        else, nORI_load = nORI; nSF_load = nSF;
        end

        fprintf('Loading and compiling bootstrapped data from BATCH')
        % Loop through batches
        for iBatch = 1:nBatches

            fprintf(' %d/%d ', iBatch, nBatches)

            iB_start = (iBatch-1)*nB_perBatch+1;
            iB_end = iBatch*nB_perBatch;

            % Define the directory to load each batch
            nameFile_boot_perBatch = sprintf('%s/%s%d/%s_batch%d_B%d_L%d%d_N%d_%d_%d%s_ORI%d_SF%d.mat', ...
                nameFolder_Data_OOD, subjName, nblocks, subjName, iBatch, nB_perBatch, iLocComb_all, flag_standEnergy, nORI, nSF, text_cut, ifamily_perF);
            load(nameFile_boot_perBatch)

            % Store data for each batch and subject
            similarityAfterMirroring_allSubj(isubj, iB_start:iB_end, :, :) = similarityAfterMirroring_allB;
            %             metrics_allSubj(isubj, iB_start:iB_end, :, :) = metrics_allB(:, :, 1:nmetrics);
            kernels2D_allSubj(isubj, iB_start:iB_end, :, :, :, :) = kernels2D_allB;
            sep_allSubj(isubj, iB_start:iB_end, :, :) = sep_allB;

            margORI_allSubj(isubj, iB_start:iB_end, :, :, :) = margORI_allB;
            margPredORI_allSubj(isubj, iB_start:iB_end, :, :, :) = margPredORI_allB;
            margParams_ORI_allSubj(isubj, iB_start:iB_end, :, :, :) = margParamsORI_allB;
            margTuningC_ORI_allSubj(isubj, iB_start:iB_end, :, :, :) = margTuningC_ORI_allB; % params2: tuning characteristics (tuningC)
            margR2ORI_allSubj(isubj, iB_start:iB_end, :, :) = margR2ORI_allB;

            margSF_allSubj(isubj, iB_start:iB_end, :, :, :) = margSF_allB;
            margPredSF_allSubj(isubj, iB_start:iB_end, :, :, :) = margPredSF_allB;
            margParams_SF_allSubj(isubj, iB_start:iB_end, :, :, :) = margParamsSF_allB;
            margTuningC_SF_allSubj(isubj, iB_start:iB_end, :, :, :) = margTuningC_SF_allB(:, :, :, 1:nTuningC_SF); % params2: tuning characteristics (tuningC)
            margR2SF_allSubj(isubj, iB_start:iB_end, :, :) = margR2SF_allB;

            if flag_ABSprefORI
                margTuningC_ORI_allSubj(:, :, :, :, 1) = abs(margTuningC_ORI_allSubj(:, :, :, :, 1));
            end

        end % ibatch
        fprintf('\n')
    end % end of isubj

    % save GROUP data (still needed for basic data, like RT)
    fprintf('\n\nSaving group data...')
    clear markers_allSubj itype
    save(nameFile_bootGroup, '*_allSubj')
    fprintf('DONE\n')
else
    fprintf('\n\nLoading group data...')
    load(nameFile_bootGroup)
    fprintf('DONE\n')
end

% plotAll_part2
