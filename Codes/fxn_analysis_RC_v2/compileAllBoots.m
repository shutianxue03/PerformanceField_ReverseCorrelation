% compile the bootstrapped outputs
clear all
clc
nB_perP = 50; % the number of bootstrappings per patch
nbatch = 20;
nB = nB_perP * nbatch;
subjList = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA', 'CS', 'DT'};
nblocks_allSubj = [200, 240, 210, 240, 220, 220, 215, 220, 205, 210, 205, 205];
subjList = { 'SX'};
nblocks_allSubj = [210];


nsubj = length(subjList);
SX_RC1_setting

for isubj = 1:nsubj
    subjName = subjList{isubj};
    nblocks = nblocks_allSubj(isubj);
    fprintf('%s (%d/%d)...', subjName, isubj, nsubj)
    
    nameIDVDdata_full = sprintf('Data_OOD/%s%d_B%d_%d_%d.mat', subjName, nblocks, nB, nfiltersOri, nfiltersSF);
    dirIDVDdata_full = dir(nameIDVDdata_full);
    
    if isempty(dirIDVDdata_full)
        SX_RC1_setting
        
        % empty containers
        dataMode = 6; SX_RC8_collectData
        
        for ibatch = 1:nbatch
            iibatch = ibatch;
            iB_start = (iibatch-1) * nB_perP+1;
            iB_end = iibatch * nB_perP;
            
            % load data
            nameIDVDdata_ = sprintf('Data_OOD/%s%d_P%d_B%d_%d_%d.mat', subjName, nblocks, ibatch, nB_perP, nfiltersOri, nfiltersSF);
            dirIDVD = dir(nameIDVDdata_);
            %         if isempty(dirIDVD), break; end
            load(nameIDVDdata_)
            
            % compile
            dataMode = 5; SX_RC8_collectData
        end
        
        % save data
        save(nameIDVDdata_full, '*_allBB', 'fitMode')
        fprintf(' DONE\n')
        % delete the patched files
        delete(sprintf('Data_OOD/%s%d_P*', subjName, nblocks))
    end
    
end % end of isubj