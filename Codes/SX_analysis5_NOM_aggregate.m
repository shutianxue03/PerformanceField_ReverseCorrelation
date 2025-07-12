
% This Script compile the model predictions and metrics derived from the
% "aggregate model" – the external variability was  the SD of all stim energy across all trials 
% (But we ended up abandoning this method according to Mike Landy's suggestion, as a model
% that predicts trial-wise response would be more mechanistic and
% informative; see NOM_trialWise)

clc, close all, warning off, format compact

addpath(genpath('Data_NOM'))
addpath(genpath('Data_compile/NOM'))
addpath(genpath('fxn_exp'))
addpath(genpath('fxn_NOM'))
addpath(genpath('fxn_analysis_RC_v2'))
addpath(genpath('XueCarrasco_JN/code'))

%% set parameters
ni = 1e3; % 
indModelsA = [1,4];
indModelsA = 1:5;
nModelsA = length(indModelsA); % should be 4, the model variation depending on randomness
templateType = 3; % template: 1=raw, 2=reconstructed, 3=mirrored, 4=mirrored&reconstructed
flagConstrainThresh = 0;% in 1 0
IVType = 1;% in 1 2 3
iModelB = 1;

% no need to vary
flagDeleteIDVDFiles = 0; % 1=delete saved files; be cautions!!
convolveType = 1; %1 2
minResample = 1;
nModelsB = 1;
flag_PatchMode = 2;
flag_standEnergy = 1;
if flag_PatchMode == 1, nameEnergySource = 'T'; else, nameEnergySource = sprintf('N%d', flag_standEnergy); end
CI_ratio = .68;

%%
nORI =29; nSF=29;
%-------------%
SX_RC1_setting
%-------------%

subjList = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA', 'CS', 'DT', 'DU', 'RC', 'SR'}; % n=15
% subjList = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD',         'HL', 'FH', 'HA',  'DT', 'CS', 'DU', 'RC', 'SR'};% n=14, no AS
% subjList = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA',  'DT', 'CS', 'DU',        'SR'};% n=14, no RC
subjList = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD',         'HL', 'FH', 'HA',  'DT', 'CS', 'DU',      'SR'};% n=13, no AS or RC
nsubj = length(subjList);

fprintf('\nnsubj=%d, nORI=%d nSF=%d\n', nsubj, nORI, nSF)

%% compile loc pairs and save all subj
% iLocComb_all = [1,8];% enter loc to compare
iLocComb_all = [6,7];
% iLocComb_all = [5,3];
iLocComb_all = [6, 7, 5, 3];
iLocComb_all = [6, 5, 3];
iLocComb_all = [6, 7, 5, 3];
nLoc = length(iLocComb_all);
if nLoc == 2, nameFileLoc = sprintf('L%d%d', iLocComb_all(1), iLocComb_all(2)); else, nameFileLoc = sprintf('Lall%d', nLoc); end

% nameFolderFig = sprintf('XueCarrasco_JN/fig/NOM/ORI%dSF%d', nORI, nSF);
nameFolderFig = sprintf('XueCarrasco_JN/fig/NOM/ORI%dSF%d/n%d', nORI, nSF, ni);
nameFolderCompile = sprintf('Data_compile/NOM/ORI%dSF%d', nORI, nSF);
if isempty(dir(nameFolderFig)), mkdir(nameFolderFig), end
if isempty(dir(nameFolderCompile)), mkdir(nameFolderCompile), end

nameFileNOM = sprintf('%s/n%d_n%d_%s_N_temp%d_B%d_nA%d_conv%d_IV%d.mat', ...
    nameFolderCompile, nsubj, ni, nameFileLoc, templateType, iModelB, nModelsA, convolveType, IVType);

% nameFileNOM = sprintf('%s/n%d_n%d_%s_B%d_nA%d.mat', nameFolderCompile, nsubj, ni, nameFileLoc, iModelB, nModelsA);

if nsubj==14 % customize file name based on which subj has been taken out
    if ~any(strcmp(subjList, 'AS')), str_noXX = 'noAS';end
    if ~any(strcmp(subjList, 'RC')), str_noXX = 'noRC';end
    nameFileNOM = sprintf('%s/n%d_%s_n%d_%s_N_temp%d_B%d_nA%d_conv%d_IV%d.mat', ...
        nameFolderCompile, nsubj, str_noXX, ni, nameFileLoc, templateType, iModelB, nModelsA, convolveType, IVType);
end
% if nModelsB>1
%     nameFileNOM = sprintf('%s/n%d_n%d_%s_N_temp%d_nA%dnB%d_conv%d_IV%d.mat', ...
%         nameFolderCompile, nsubj, ni, nameFileLoc, templateType, nModelsA, nModelsB, convolveType, IVType);
% end

dirModel = dir(nameFileNOM);

if isempty(dirModel)
    
    NOM_data_metrics_allSubj = nan(nsubj, nLoc, nModelsA, nModelsB, ni, 8); % 8 metrics [data from the test group]
    NOM_pred_metrics_allSubj  = NOM_data_metrics_allSubj; % 8 metrics [pred]
    NOM_params_est_allSubj = nan(nsubj, nLoc, nModelsA, nModelsB, ni, 3); % estimated params, hard code 3
    NOM_nLL_allSubj = nan(nsubj, nLoc, nModelsA, nModelsB, ni); % estimated params
    NOM_h_allSubj = nan(nsubj, nLoc, nModelsA, nModelsB, ni, 4,2); % what is this?
    NOM_maxChannel_allSubj = nan(nsubj, nLoc, nModelsA, nModelsB, 2,2);
    NOM_IC_allSubj = nan(nsubj, nLoc, nModelsA, nModelsB, ni, nIC);
    
    for iiLoc = 1:nLoc
        iLocComb = iLocComb_all(iiLoc);
        fprintf('\n\n%s [L%d]...', namesLocComb{iLocComb}, iLocComb)
        
        for isubj = 1:nsubj
            subjName = subjList{isubj};
            fprintf('\n  %s %d/%d...\n', subjName, isubj, nsubj)
            
            nameFolderNOM = sprintf('Data_NOM/ORI%dSF%d/%s/L%d', nORI, nSF, subjName, iLocComb);
            
            for iModelA = indModelsA%1:nModelsA
%                 for iModelB = 2%1:nModelsB
                    if iModelB>1, nNoise=2; else, nNoise=3; end
                    
                    nameFileModelIDVD = sprintf('%s/n%d_N_A%dB%d_temp%d_min%d_thresh%d_conv%d_IV%d.mat', ...
                        nameFolderNOM, ni, iModelA, iModelB, templateType, minResample, flagConstrainThresh, convolveType, IVType);
%                     nameFileModelIDVD = sprintf('%s/n%d_A%dB%d.mat', nameFolderNOM, ni, iModelA, iModelB);
                    dirFile = dir(nameFileModelIDVD);
                    
                    fprintf('    A%dB%d', iModelA, iModelB)
                    
                    if isempty(dirFile)
                        fprintf(' [Not exist]\n')
                    else
                        load(nameFileModelIDVD)
                        fprintf(' Loaded\n')
                    
                    % compile data
                    NOM_data_metrics_allSubj(isubj, iiLoc, iModelA, iModelB, :, :) = data_metrics_allB;
                    NOM_pred_metrics_allSubj(isubj, iiLoc, iModelA, iModelB, :, :) = pred_metrics_allB;
                    NOM_params_est_allSubj(isubj, iiLoc, iModelA, iModelB, :, 1:nNoise) = params_est_allB;
                    NOM_nLL_allSubj(isubj, iiLoc, iModelA, iModelB, :) = nLL_allB;
                    NOM_h_allSubj(isubj, iiLoc, iModelA, iModelB, :, :, :) = h_allB; % 1=passed swtest for 4 transforms: {'Raw'}    {'Log'}    {'Sqrt'}    {'boxcox'}
                    NOM_IC_allSubj(isubj, iiLoc, iModelA, iModelB, :, :) = IC_allB; % nsubj x nLoc x nModelsA x nModelsB x ni x nIC
                    end
                    % intermediate outputs
                    %                     maxChannel_med_allB = nan(ni, 2, 2); % PRS/ABS, ORI/SF
                    %                     for ii = 1:ni
                    %                         maxChannel_med_allB(ii, 1, :) = median(data_allB{ii}.imax{1});
                    %                         maxChannel_med_allB(ii, 2, :) = median(data_allB{ii}.imax{2});
                    %                     end
                    %                     NOM_maxChannel_allSubj(isubj, iiLoc, iModelA, iModelB, :, :) = squeeze(median(maxChannel_med_allB, 1));
                    
                    if flagDeleteIDVDFiles
                        delete(nameFileModelIDVD)
                    end
%                 end % iModelB
            end % iModelA
        end % iLoc
    end % isubj
    fprintf('\nDONE\n')
    
    % save
    save(nameFileNOM, 'NOM*_allSubj')
else, load(nameFileNOM), fprintf('\n\n === Loaded === \n\n\n')
end

%% plot Fig 1-3
flag_plotCI = 0; close all
%-----------%
modelPlot0
%-----------%

%% Fig 4. compare candidate models
% if want to run this section again, you need to run the section above '%% plot Fig 1-3'
flag_plotIDVD = 0; 
iIC = 0; % 0=nLL; 1=AIC; 2=AICc; 3=BIC % the resulting dalta are the same
nPairs = 2; % for bonferroni correction
str_tail = 'left'; % assume A<B
close all
%-----------%
modelPlot4_withB4
%-----------%

%-----------%
% modelPlot4
%-----------%

%% Fig 5. prop. of tests passing the normality test (swtest)
% close all
% for iLocComb = 1:nLoc
% modelPlot5
% end

%% Fig 6. imax of the channel
% at which ORI/SF the max of the dotProduct is found
% should be consistent with the peak of the tuning function
% close all
% for iiLoc = 1:nLoc
%     iLocComb = iLoc_all(iiLoc);
%     modelPlot6
% end
