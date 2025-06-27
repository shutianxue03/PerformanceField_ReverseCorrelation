
% For each observer, this script loads their data (processed on OOD), plot
% idvd data and for all observers, plot group data.
% path(pathdef)
clc
close all
clear all
warning off
format compact

% addpath(genpath('XueCarrasco')) % this takes a little while
addpath(genpath('fxn_MC'))
addpath(genpath('Data_OOD'))
addpath(genpath('Data_compile/RC'))
addpath(genpath('fxn_exp'))
addpath(genpath('fxn_analysis_RC_v2'))
addpath(genpath('fxn_simulation'))
addpath(genpath('fxn_RCplot'))
addpath(genpath('XueCarrasco_JN/code'))

%% basic settings
subjList = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA',  'DT', 'CS', 'DU', 'SR', 'RC'};
nblocks_allSubj = [200, 240, 210, 240, 220, 220, 215, 220, 205, 210, 205, 205, 205, 195, 205];
% subjList = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD',       'HL', 'FH', 'HA',  'DT', 'CS', 'DU',      'SR'};
% nblocks_allSubj = [200, 240, 210, 240, 220, 220,       220, 205, 210, 205, 205, 205,       195];

% subjList = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA',  'DT', 'CS', 'DU', 'SR'};
% nblocks_allSubj = [200, 220, 210, 220, 220, 220, 215, 220, 205, 210, 205, 205, 205, 195];

nsubj = length(subjList);

global flag_cutMapping ifamily_perF

nModelsA=5;

%---------------%
SX_RC1_setting
%---------------%
flag_ABSprefORI = 1;%input('     >>> take abs of pref ORI? (1=YES, 0=NO): '); % 1=take abs, meanig distance from the signal; 0=the sign (+ or -) indicates bias
flag_cutMapping = 0; % 1=load original data (ORI29 SF29) and cut the edge
flag_mirrorMapping = 1; % 1=mirror 2D mapping over 0 deg ori
nB_perBatch = 1e2;
nBatches = 10;
nB = nB_perBatch * nBatches;
ni = 1e3; 
flag_PatchMode = 2; % energy derives from target; 2=noise patch
flag_block200 = 0; % 1; all observers are forced to have 200 blocks; 0=no
flag_standEnergy = 1; if flag_PatchMode==1, flag_standEnergy =1; end

% determine the list of subjects
nsubj = length(subjList);
namesEnergySource = {'TARGET', 'NOISE'};
if flag_PatchMode == 1, namePatchMode = 'T';
else, namePatchMode = sprintf('N%d', flag_standEnergy);
end


%% define tuning fxn fitting params
edit plotAll_part2.m
edit plotAll_part1.m
%---------------%
SX_RC1_setting
%---------------%
fprintf('\n\nnORI=%d nSF=%d\n\n', nORI, nSF)

iLocComb_all_all = {[6,7], [5,3]};
% [1,8]: fovea vs. periphery
% [2,4]: LHM vs. RHM (also plot D vs. nonD) % run this before [6,7]!!
% [6,7]: HM vs. VM (also plot HM (D) vs. VM)
% [5,3]: LVM vs. UVM

ill = input('     >>> Which location to compare (1=L18, 2=L67, 3=L53, 4=L24, 5=L47): ');
iLocComb_all = iLocComb_all_all{ill};% delete

ifamily_perF = input('     >>> What is the model family for ORI and SF (e.g., [8,2]): '); % ORI: Gaussian with different mean; SF: double peak

nTuningC_ORI = length(namesTunC_unit_perF{ifamily_perF(1), 2}); % preferred ori, peak amp, bandwidth, bottom, // pred_half
nTuningC_SF = length(namesTunC_unit_perF{ifamily_perF(2), 2}); % peakSF, peak amp., bandwidth_full (octave), baseline // bandwidth_L, bandwidth_R, bottom, pred_half, full bandwidth (linear scale)

paramInd_perF{1} = ones(1,length(namesParams_all{ifamily_perF(1)})); %
paramInd_perF{2} = ones(1,length(namesParams_all{ifamily_perF(2)}));

fprintf('\nnsubj = %d, nB = %d\n\n%s vs. %s\nEnergy derived from %s patch (standardization=%d)\n\nORI: M%d %s [%s]\nSF: M%d %s [%s]\n', ...
    nsubj, nB, ...
    namesLocComb{iLocComb_all(1)}, namesLocComb{iLocComb_all(2)}, ...
    namesEnergySource{flag_PatchMode}, flag_standEnergy, ...
    ifamily_perF(1), namesFamily_all{ifamily_perF(1)}, num2str(paramInd_perF{1}), ...
    ifamily_perF(2), namesFamily_all{ifamily_perF(2)}, num2str(paramInd_perF{2}))

%% group file name
if flag_mirrorMapping, text_cut = '_m'; else, text_cut = ''; end

nameFolderCompile = 'Data_compile/RC';
nameFileBootGroup = sprintf('%s/%d%d%s/n%d_B%d_L%d%d_%d_%d%s_ORI%dSF%d.mat', ...
    nameFolderCompile, nORI, nSF, text_cut, nsubj, nB, iLocComb_all, nORI, nSF, text_cut, ifamily_perF(1), ifamily_perF(2));

dirBootGroup = dir(nameFileBootGroup);

if isempty(dirBootGroup)
    
    % empty holders
    %     if nmetrics == 9, nmetrics=8; end
    %     metrics_allSubj = nan(nsubj, nB,  nLoc2, nmetrics);
    %     cst_allSubj = nan(nsubj, nB,  nLoc2);
    %     RT_allSubj = nan(nsubj, nB,  nLoc2);
    pYES_tgt_allSubj = nan(nsubj, nB,  nLoc2, ntypes, nbins_e);
    ebin_tgt_allSubj = nan(nsubj, nB,  nLoc2, ntypes, nbins_e);
    similarityAfterMirroring_allSubj = nan(nsubj, nB,  ntypes, nLoc2);
    kernels2D_allSubj = nan(nsubj, nB,  nLoc2, ntypes, nORI, nSF);
    sep_allSubj = nan(nsubj, nB,  nLoc2, ntypes);
    
    margORI_allSubj = nan(nsubj, nB,  nLoc2, ntypes, nORI);
    margORIraw_allSubj = nan(nsubj, nB,  nLoc2, ntypes, nORI);
    margPredORI_allSubj = nan(nsubj, nB,  nLoc2, ntypes, nORI);
    margParams_ORI_allSubj = nan(nsubj, nB,  nLoc2, ntypes, length(paramInd_perF{1}));
    margTuningC_ORI_allSubj = nan(nsubj, nB,  nLoc2, ntypes, nTuningC_ORI);
    margR2ORI_allSubj = nan(nsubj, nB,  nLoc2, ntypes);
    
    margSF_allSubj = nan(nsubj, nB,  nLoc2, ntypes, nSF);
    margPredSF_allSubj = nan(nsubj, nB,  nLoc2, ntypes, nSF);
    margParams_SF_allSubj = nan(nsubj, nB,  nLoc2, ntypes, length(paramInd_perF{2}));
    margTuningC_SF_allSubj = nan(nsubj, nB,  nLoc2, ntypes, nTuningC_SF);
    margR2SF_allSubj = nan(nsubj, nB,  nLoc2, ntypes);
    
    %
    for isubj = 1:nsubj
        
        subjName = subjList{isubj};
        nblocks = nblocks_allSubj(isubj);
        fprintf('\n%s (%d/%d)...\n', subjName, isubj, nsubj)
        
        % load behav
        %         nameBehavMeas = sprintf('Data_OOD/%s%d/%s_behavMeas.mat', subjName, nblocks, subjName);
        %         fprintf('Loading behav data...')
        %         load(nameBehavMeas)
        %         fprintf(' DONE\n')
        
        % load boot data files
        if flag_cutMapping, nORI_load = nan; nSF_load = nan;
        else, nORI_load = nORI; nSF_load = nSF;
        end
        
        %         nameFileBoot = sprintf('Data_OOD/%s%d/%s_B%d_L%d%d_N%d_%d_%d%s_ORI%d_SF%d.mat', ...
        %             subjName, nblocks, subjName, nB, iLocComb_all, flag_standEnergy, nORI, nSF, text_cut, ifamily_perF);
        
        fprintf('Loading and compiling BOOTs...\n')
        %         load(nameFileBoot)
        
        for ibatch = 1:nBatches
            iB_start = (ibatch-1)*nB_perBatch+1;
            iB_end = ibatch*nB_perBatch;
            nameFileBoot_perBatch = sprintf('Data_OOD/%s%d/%s_batch%d_B%d_L%d%d_N%d_%d_%d%s_ORI%d_SF%d.mat', ...
                subjName, nblocks, subjName, ibatch, nB_perBatch, iLocComb_all, flag_standEnergy, nORI, nSF, text_cut, ifamily_perF);
            
            fprintf(' #%d ', ibatch)
            %             if isempty(dir(nameFileBoot_perBatch)), fprintf('NOT exist\n'), continue, else, fprintf('\n'), end
            load(nameFileBoot_perBatch)
            
            % compile all subj
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
            
            fprintf('.')
        end % ibatch
    end % end of isubj
    
    % save GROUP data (still needed for basic data, like RT)
    fprintf('\n\nSaving...')
    clear markers_allSubj itype
    save(nameFileBootGroup, '*_allSubj')
    fprintf('DONE\n')
else
    fprintf('Loading...')
    load(nameFileBootGroup)
    fprintf('DONE\n')
end

% plotAll_part2
