% compute evaluated data across all observers
clc
close all

addpath(genpath('fxn_exp/'))
addpath(genpath('fxn_analysis_titration/'))
addpath(genpath('fxn_analysis_RC/'))
addpath(genpath('palamedes/'))
addpath(genpath('Data/'))

%%
% subjList = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA', 'CS', 'DT'};
subjList = {'QS'}; 
nsubj = length(subjList);

draw_mode = 1;
exp_mode = 0;
data_allSubj = nan(5,11,nsubj); % 5 loci x 11 measurements x nsubj
ratioSelected = .7; % proportion of trials selected for computing accuracy
fprintf('Last %d%% trials are evaluated.\n', round(ratioSelected*100))

locNames = {'center', 'left', 'upper', 'right', 'lower'};
nLoc = length(locNames);

% thresh_final_allSubj =  [ 0.36, 0.5, 0.715, 0.5, 0.55];

thresh_final_allSubj = [
    [ 0.315, 0.36, 0.37, 0.365, 0.35];            % 1. SX
     [.34 .49 .505 .42 .41];                            % 2. SP
    [ 0.375, 0.43, 0.445, 0.41, 0.345];          % 3. LS
     [ 0.375, 0.4325, 0.405, 0.425, 0.392];  % 4. RE
      [ 0.365, 0.455, 0.715, 0.455, 0.545];     % 5. YK
      [ 0.35, 0.475, 0.5, 0.53, 0.515];             % 6. HL
      [ 0.35, 0.345, 0.39, 0.365, 0.425]];      % 7. DT

thresh_final_allSubj = [
    [ 0.2, 0.235, 0.26, 0.235, 0.245]; % SX
    [ 0.34, 0.43, 0.5, 0.43, 0.46]; % FH
    [ 0.365, 0.465, 0.6, .47, 0.49]; % YK
    [ 0.36, 0.42, 0.43, 0.42, 0.425]; % ME (Marissa)
    [ 0.38 0.395, 0.415, 0.395, 0.385]; % MD 
    [ 0.33, 0.34, 0.39, 0.34, 0.35];%LL
    [ 0.4573, 0.4507, 0.3593, 0.3666, 0.4943]; % CS (echo)
];clear all


%% plot individual data
for isubj = 1:nsubj
    subjName = subjList{isubj};
    folderName = ['Data/', subjName, '/'];
    plotErrorBar = 0;
    fprintf('%s (%d/%d) in progress ...', subjName, isubj, nsubj)
    
    % 1. load data
    dataFile_dir = dir([folderName,subjName,'_stair_','*']);
    if isempty(dataFile_dir), error('ALERT: No titration data of this subject.'),
    else, load(dataFile_dir(end).name,'record'), fprintf(' titration data loaded.')
    end
    
    % 2. plot data
    drawEachSubj
    set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
    
    % 3. collect data
    data_allSubj(:,:,isubj) = [thresh_log, thresh_log_mean, acc_prs', acc_abs', acc_all', ...
        RT_prs', RT_abs', RT_all', dprime', criterion'];
    
    clear thresh_final_log
end

if nsubj > 1, save(sprintf('Data/titrationMat_n%d', nsubj), 'data_allSubj'), end

%% plot average data 
% if nsubj > 1, drawAllSubj, end




