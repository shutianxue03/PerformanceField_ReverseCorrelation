clc
close all

addpath(genpath('Data_OOD'))
addpath(genpath('fxn_model'))
addpath(genpath('fxn_simulation'))
addpath(genpath('fxn_RCplot'))
clear data_allSubj pred_med_allSubj

%% SETTING
subjList = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA', 'CS', 'DT'};
nblocks_allSubj = [200, 240, 210, 240, 220, 220, 215, 220, 205, 210, 205, 205];

% subjList = {'YK', 'SP', 'SX', 'LS', 'MD', 'AS', 'HL', 'FH', 'HA', 'CS', 'DT'};
% nblocks_allSubj = [200, 240, 210, 240, 220, 215, 220, 205, 210, 205, 205];


% subjList = {'YK', 'SX'};
% nblocks_allSubj = [200, 210];
% 
% subjList = {'SX'};
% nblocks_allSubj = 210;

modelVersion = 3; % find details of this version in the var 'readme'
renameMode = 0; % 
iLoc_all = [1,8];
nparams_model = 3; 

%%
nsubj = length(subjList);
nLoc = length(iLoc_all);

%% plot and publish
publishOptions = struct('format','pdf','outputDir','publishedPDFs/', 'showCode', boolean(0));
if nsubj == 1
    subjName = subjList{1};
    publishName = sprintf('publishedPDFs/%s_model%d', subjName, modelVersion);
else
    publishName = sprintf('publishedPDFs/n%d_model%d', nsubj, modelVersion);
end

publish('modelPlot_all', publishOptions);
close all
movefile('publishedPDFs/modelPlot_all.pdf', publishName);
disp('PDF published.')



