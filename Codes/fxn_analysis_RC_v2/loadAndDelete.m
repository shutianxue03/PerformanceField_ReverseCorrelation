clc
clear all

isubj = 12;

subjList = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA', 'CS', 'DT'};
nblocks_allSubj = [200, 240, 210, 240, 220, 220, 215, 220, 205, 210, 205, 205];
subjName =  subjList{isubj};
nblocks = nblocks_allSubj(isubj);
% load
load(sprintf('Data_OOD/%s%d/%s_behavMeas.mat', subjName, nblocks, subjName))
% clear
clear isubj subjList nblocks_allSubj nsubj
% save
save(sprintf('Data_OOD/%s%d/%s_behavMeas.mat', subjName, nblocks, subjName))
clear all
fprintf('DONE')
