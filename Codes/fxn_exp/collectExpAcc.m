
% subjName = input('enter subj name: ');
subjList = {'SX', 'AD', 'NH', 'ZZ'};
nsubj = length(subjList);
for isubj = 1:nsubj
    subjName = subjList{isubj};
    folderName = ['Data/', subjName, '/'];
    load([folderName, subjName, 'params1.mat'])
    delete([folderName, subjName, '_ThreshAccRecord.mat'])
    
    threshT_allSess = nan(1, 5);
    accT_allSess = nan(1, 5);
    accExp_allSess = nan(1, 5);
    
    %% save titration
    stairFileDir = dir([folderName, subjName, '_stair*']); % find out all possible files
    load(stairFileDir(end).name)
    saveThreshAcc
    
    %% save exp accuracy
    expFileDir = dir([folderName, subjName, '_exp*']); % find out all possible files
    if isempty(expFileDir), error('ALERT: no exp data for %s !!!', subjName), end
    nblocks = length(expFileDir); % number of blocks finished
    for ib = 4:4:nblocks
        load(expFileDir(ib).name)
        saveThreshAcc
    end
end