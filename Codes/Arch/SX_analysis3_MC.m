
% Last updated on 07/14/2025 by Shutian Xue 

% This script performs model comparison (MC) analysis for the
% psychophysical data collected in the RC lab. It compares different
% families of models (e.g., Gaussian, Difference of Gaussians, Log
% Parabola, Truncated Log Parabola) based on their fit to the data
% across different observers and conditions. The analysis includes
% cross-validation, leave-one-out cross-validation, and information
% criteria (AIC, AICc, BIC) to determine the best model for each
% feature (Orientation, Spatial Frequency) and family of models.


%%%%%%%%%
%  Inputs
%%%%%%%%%
% ifeature: 1=ORI, 2=SF
% ifamily_all: ORI: 1 (Gaussian) or 8 (DoG); SF: 2 (log parabola), 3 (truncated log parabola)
% MCmode: 1=10 fold CV, 2=LOOCV, 3=Info criterion
% nIC: when MCmode = 3, 1=AIC, 2=AICc, 3=BIC

%%%%%%%%%
% Outputs
%%%%%%%%%
% For each ifeature, ifamily and MCmode:
% Fig 1: [top] dev of each model averaged across subj
%           [bottom] freq of the model to be chosen as the best model
% Fig 2: data + pred for each observer (two fitting lines: the best idvd model and second, the best group model)

% For each ifeature and MCmode:
% Fig 3: bar plots of dev of each ifamily for each MCmode

clear all, close all, clc
addpath(genpath('staircasecode_MJ'))
addpath(genpath('fxn_analysis_RC_v2'))
addpath(genpath('fxn_MC'))
addpath(genpath('fxn_exp'))
% addpath(genpath('XueCarrasco_JN/code'))

datetime('now')

%%
clc
ifamily_all_perF = {[1,8], [2, 3]}; % models to compare
iLocPairs_all = [1,8; 6,7; 5,3]; %location pairs
nLocPairs = size(iLocPairs_all, 1);

% Define flags
flag_cutMapping = 0; %input('       >>> Whether cut the edges (1=CUT, 0=NO): ');
flag_mirrorMapping=1; %input('       >>> Whether mirror the kernel mapping (1=mirror, 0=NO): ');
MCmode = 1; %input('       >>> What model comparison mode (1=CV, 2=LOOCV, 3=iC): ');
flag_plotMC = 1; %input('       >>> Make plots? (1=YES, 0=NO): ');
flag_devErrorbar = 0;%input('       >>> Plot errorbars in dev plots? (1=YES, 0=NO): ');
flag_plotMC_IDVD  = 0;%input('       >>> Plot idvd data in family comparison plots? (1=YES, 0=NO): ');

% nORI = 29; nSF = 29;
%-------------%
SX_RC1_setting
%-------------%
subjList = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA',  'DT', 'CS', 'DU', 'SR', 'RC'};
nblocks_allSubj = [200, 240, 210, 240, 220, 220, 215, 220, 205, 210, 205, 205, 205, 195, 205];

% subjList =            {'YK', 'SP', 'SX', 'LS', 'RE', 'MD',         'HL', 'FH', 'HA',  'DT', 'CS', 'DU',        'SR'}; % no AS or RC, n=13
% nblocks_allSubj = [200, 240, 210, 240, 220, 220,        220, 205, 210, 205, 205, 205,      195];
nsubj = length(subjList);

flag_PatchMode = 2; % energy derives from target; 2=noise patch
if flag_PatchMode == 1, namePatchMode = 'T'; else, namePatchMode = 'N'; end
iType = 2;
close all

flag_plotIDVDdev = 0; % 1=plot idvd dev as grey line in Fig 1 (MCplot1)
flag_plotComp = 1; % 1=plot comparison of dev/IC
flagPlotPred = 1; % 1=plot idvd prediction, 0=do NOT plot
flag_fancyXticks = 0; %1=plot candidate index (1=free to vary; 0= shared)

%%
for iiLocPair = 1:nLocPairs
    
    iLocComb_all = iLocPairs_all(iiLocPair, :);
    %     dev_bestGroup_perModePerFamily = cell(8, 3, nIC); % 8 = number of families, 3 = number of MCmodes, 3=nIC
    if flag_mirrorMapping, nameM = '_m'; else, nameM = ''; end
    if flag_cutMapping, nameC = 'c'; else, nameC = ''; end
    
    for iFeature = 1:2
        ifamily_all = ifamily_all_perF{iFeature};
        
        % Define Figure folder
        nameFolder_Fig_MC = sprintf('%s/MC/%d%d%s/L%d%d/%s/', nameFolder_Figures, nORI, nSF, nameM, iLocComb_all, namesFeature{iFeature});
        nameDir_Fig_MC = dir(nameFolder_Fig_MC);
        if isempty(nameDir_Fig_MC), mkdir(nameFolder_Fig_MC), end
        
        ibestGroup_perModePerFamily = cell(1,length(ifamily_all));
        
        % Loop through model families
        for ifamily = ifamily_all
            %% compile idvd file
            nameFile_MCallSubj = sprintf('%s/Data_MC/%d%d%s/n%d_L%d%d_Family%d_mode%d.mat', nameFolder_Data, nORI, nSF, nameM, nsubj, iLocComb_all, ifamily, MCmode);
            
            %% model setting
            namesParams = namesParams_all{ifamily};
            nparams_full = length(namesParams);
            nCands = 2^nparams_full;
            paramInd_all = fxn_getParamInd(nparams_full);
            
            if MCmode == 3, nIC_ = nIC; else, nIC_ = 1; end
            
            if isempty(dir(nameFile_MCallSubj))
                %% Preallocate variables
                irank_allSubj = nan(nsubj, nIC_, nCands); % the sorted model index (best to worst model);
                iBest_allSubj = nan(nsubj, nIC_);
                dev_allSubj = nan(nsubj, nCands, 3);
                params_est_allSubj = cell(nsubj, nCands);
                
                for isubj = 1:nsubj
                    fprintf('%s (%d/%d)', subjList{isubj}, isubj, nsubj)
                    
                    % Data folder
                    nameFolder_MC = sprintf('%s/Data_MC/%d%d%s%s/%s/L%d%d/', nameFolder_Data, nORI, nSF, nameM, nameC, subjList{isubj}, iLocComb_all);
                    nameFile_MC = sprintf('%sFamily%d_mode%d.mat', nameFolder_MC, ifamily, MCmode);
                    if isempty(dir(nameFile_MC)), fprintf('\n********************\n NOT EXIST: %s\n********************\n', nameFile_MC)
                    else, load(nameFile_MC), fprintf('... Loaded\n')
                    end
                    
                    for iIC = 1:nICs
                        irank_allSubj(isubj, :, :) = irank_allCand; % the sorted model index (best to worst model);
                        iBest_allSubj(isubj, :) = iBest_allCand;
                    end
                    dev_allSubj(isubj, :, :) = dev_allCand;
                    params_est_allSubj{isubj} =params_est_allCand;
                    
                end % isubj
                
                save(nameFile_MCallSubj, 'ifeature', 'ifamily',  'MCmode', 'paramInd_all', '*_allSubj')
            else, load(nameFile_MCallSubj), fprintf('Loaded\n')
            end
            
            %% Fig 1-2
            %-------------------%
            MCplot0_allMode
            %-------------------%
            close all
            ibestGroup_perModePerFamily{ifamily} = num2str(paramInd_all(iBest_group, :));
            
        end % ifamily
        
        % Fig 3. compare across families
        if length(ifamily_all)>1
            %-------------------%
            MCplot3_compFamily
            %-------------------%
        end
    end % ifeature
    
end % iLoc_all_all

close all

