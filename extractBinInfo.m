% 1. extract the number of days each observer took, and the days between sessions
% 2. decide how to bin data
% 3. based on the binned data, calculate (1) the number of days each bin
% took (2) the number of days between bins

clc, close all, warning off, format compact

addpath(genpath('Data'))
addpath(genpath('fxn_analysis_RC_v2'))
addpath(genpath('XueCarrasco_JN/code'))

%%
clc,close all
% subjList_full              = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA',  'DT', 'CS', 'DU', 'SR', 'RC'};
% nblocks_allSubj_full = [200, 240, 210, 240, 220, 220, 215, 220, 205, 210, 205, 205, 205, 195, 205];
% indSubj = 1:15;
% indSubj = [1:6, 8:14]; % no AS or RC
% indSubj = 1:11;
% subjList = subjList(indSubj);
% nblocks_allSubj = nblocks_allSubj(indSubj);

subjList =             {'YK', 'SP', 'LS', 'RE', 'MD', 'HL', 'FH', 'HA',  'DT', 'CS', 'SR',    'SX', 'DU'}; % the last two are authors
nblocks_allSubj = [200, 240, 240, 220, 220, 220, 205, 210,   205, 205, 195,    210, 205];
nsubj = length(subjList);

nLoc = 5; % each session contains 5 blocks, each testing one loc once; subj always finished multiples of 5 blocks per day
nBins = 3;
binStrategy_all = {'equal', 'algorithm', 'manual_SX', 'manual_DT'};

%%
% Extract date information
% outputs:
%       days_total_allSubj: nsubj x 1, total duration of data collection in days for each subj
%       date_perSess_allSubj: nsubj x 1 / nSess x 1, the date (str, in datetime format) on which each session was completed
%       days_diffFromLastSess_allSubj: nsubj x 1 / [nSess-1 x 1], the gap from last session in days

days_total_allSubj = nan(nsubj, 1);
date_perSess_allSubj = cell(nsubj, 1);
gap_diffFromLastSess_allSubj = cell(nsubj, 1);

for isubj = 1:nsubj
    
    subjName = subjList{isubj};
    nblocks = nblocks_allSubj(isubj);
    filesDir = dir(sprintf('Data/%s/%s_exp_B*', subjName, subjName));
    if strcmp(subjName, 'SP'), filesDir = filesDir(81:end); end
    assert(length(filesDir) == nblocks)
    
    nSess = nblocks/nLoc;
    
    % extract the date of each session
    date_perSess = cell(nSess, 1);
    gap_diffFromLastSess = nan(nSess-1, 1);
    for iSess = 1:nSess
        iblock = (iSess-1)*nLoc+1;
        val_year = str2double(filesDir(iblock).name(15:18));
        val_month = str2double(filesDir(iblock).name(19:20));
        val_date= str2double(filesDir(iblock).name(21:22));
        date_perSess{iSess} = datetime(val_year, val_month, val_date);
        
        if iSess>1
            gap_diffFromLastSess(iSess-1) = round(days(date_perSess{iSess} -date_perSess{iSess-1}));
        end
    end % iSess
    
    date_perSess_allSubj{isubj} = date_perSess;
    gap_diffFromLastSess_allSubj{isubj} = gap_diffFromLastSess;
    
    % Calculate the total duration (in days)
    d1 = date_perSess{1};
    dLast = date_perSess{end};
    days_total_allSubj(isubj) = round(days(dLast-d1));
    
end

%% combine algorithm-based and manual search
clc
nSess_perBin_allSubj = cell(nsubj, 1);
dur_perBin_allSubj = nSess_perBin_allSubj;
gap_btwBin_allSubj = nSess_perBin_allSubj;
iSess_start_perBin_allSubj = nSess_perBin_allSubj;
iSess_end_perBin_allSubj = nSess_perBin_allSubj;

flag_binStrategy = input('       >>> 0=equal; 1=Algorithm-generated strategy; 2=Manually decided strategy (SX); 3=manual (DT)');
binStrategy = binStrategy_all{flag_binStrategy+1};
switch flag_binStrategy
    case 0
        nSess_perBin_allSubj = {[13, 13, 13]
            [13, 13, 13]
            [13, 13, 13]
            [13, 13, 13]
            [13, 13, 13]
            [13, 13, 13]
            [13, 13, 13]
            [13, 13, 13]
            [13, 13, 13]
            [13, 13, 13]
            [13, 13, 13]
            [13, 13, 13]
            [13, 13, 13]};
        
    case 1
        nCombBest = 5;
        nSess_perBin_min = 10; % the min. number of sessions each bin contains
    case 2
        nSess_perBin_allSubj = {[12, 14, 14] % YK
            [19, 15, 14] % SP
            [16, 23, 9] % LS
            [15, 14, 15] % RE
            [15, 15, 14] % MD
            [8, 20, 16] % HL
            [12, 13, 16] % FH
            [14, 14, 14] % HA
            [10, 16, 15] % DT
            [13, 11, 17] % CS
            [13, 13, 13] % SR
            [12, 13, 17]  % SX
            [11, 9, 21] }; % DU
end

for isubj = 1:nsubj
    clc
    subjName = subjList{isubj};
    nblocks = nblocks_allSubj(isubj);
    nSess = nblocks/nLoc;
    date_perSess = date_perSess_allSubj{isubj};
    gap_diffFromLastSess = gap_diffFromLastSess_allSubj{isubj};
    nSess_perBin_suggested = round(nSess/nBins);
    
    if flag_binStrategy==1
        %%%%%% ALGORITHM-BASED %%%%%%
        % the starting session of all but the first bin
        iSess_start_allComb = combvec((nSess_perBin_min+1):nSess, (nSess_perBin_min+1)*2:nSess);
        nComb = size(iSess_start_allComb, 2);
        
        val_allComb = nan(nComb, 1);
        for iComb = 1:nComb
            iSess_start_allBins = iSess_start_allComb(:, iComb);
            if iSess_start_allBins(1)>=iSess_start_allBins(2)
                continue
            else
                %------------------------%
                [nSess_perBin, dur_perBin, gap_btwBin] = fxn_getBinInfo(nBins, nSess, date_perSess, iSess_start_allBins, gap_diffFromLastSess);
                %------------------------%
                
                if any(nSess_perBin<nSess_perBin_min)
                    continue
                else
                    val_allComb(iComb) = std(nSess_perBin) - sum(gap_btwBin) + sum(dur_perBin);
                end
            end
        end % iComb
        
        [val_min, iComb_best] = sort(val_allComb);
        
        fprintf('\n============================')
        fprintf('\n%s nSess=%d Total dur = %d days [max. gap=%d] (suggested nSess per bin): %d', ...
            subjName, nSess, days_total_allSubj(isubj), max(gap_diffFromLastSess), nSess_perBin_suggested)
        fprintf('\n============================\n')
        
        for iCombBest = 1:nCombBest
            fprintf('\n*** [#%d] Val = %.2f***\n', iCombBest, val_min(iCombBest))
            
            %%%%%%%%%%%%
            iSess_start_allBins = iSess_start_allComb(:, iComb_best(iCombBest));
            [nSess_perBin, dur_perBin, gap_btwBin] = fxn_getBinInfo(nBins, nSess, date_perSess, iSess_start_allBins, gap_diffFromLastSess);
            %%%%%%%%%%%%
            
            fprintf('Number of sessions per bin: [%s] SD=%.1f\n', num2str(nSess_perBin'), std(nSess_perBin))
            fprintf('Gaps (in days) between bins: [%s]\n', num2str(gap_btwBin'))
            fprintf('Duration (in days) per bin: [%s] \n', num2str(dur_perBin'))
        end
        iSess_start_allBins = iSess_start_allComb(:, iComb_best(1));
        %-------------------------%
        [nSess_perBin, dur_perBin, gap_btwBin] = fxn_getBinInfo(nBins, nSess, date_perSess, iSess_start_allBins, gap_diffFromLastSess);
        %-------------------------%
    else %%%%%% MANUAL %%%%%%
        nSess_perBin = nSess_perBin_allSubj{isubj}.';
        nSess = sum(nSess_perBin);
        iSess_start_allBins = cumsum(nSess_perBin);
        iSess_start_allBins = iSess_start_allBins(1:nBins-1)+1;
        %-------------------------%
        [nSess_perBin_, dur_perBin, gap_btwBin] = fxn_getBinInfo(nBins, nSess, date_perSess, iSess_start_allBins, gap_diffFromLastSess);
        %-------------------------%
        assert(sum(nSess_perBin_-nSess_perBin)==0)
    end % if flagLoadBinStrategy
    
    nSess_perBin_allSubj{isubj} = nSess_perBin;
    dur_perBin_allSubj{isubj} = dur_perBin;
    gap_btwBin_allSubj{isubj} = gap_btwBin;
    iSess_start_perBin_allSubj{isubj} = [1, cumsum(nSess_perBin(1:nBins-1)')+1];
    iSess_end_perBin_allSubj{isubj} = cumsum(nSess_perBin');
end % isubj

save(sprintf('Data_compile/BinStrategy/BinInfo_n%d_nBins%d_%s.mat', nsubj, nBins, binStrategy), '*_allSubj')

%% display the info of the best binning strategy (either by algorithm or by manual selection)
clc
% nBins = 3;
% nLoc = 5;
% subjList_full = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA',  'DT', 'CS', 'DU', 'SR', 'RC'};
% nblocks_allSubj_full = [200, 240, 210, 240, 220, 220, 215, 220, 205, 210, 205, 205, 205, 195, 205];
colors_bin = [.75, .75, .75; .5, .5, .5; 0, 0, 0];

for flag_binStrategy=0:2
    clc
    binStrategy = binStrategy_all{flag_binStrategy+1};
    if flag_binStrategy == 3, indSubj = 1:11; else, indSubj = [1:6, 8:14]; end
%     subjList = subjList_full(indSubj);
%     nblocks_allSubj = nblocks_allSubj_full(indSubj);
%     nsubj = length(subjList);
    
    nameFigFolder_bin = sprintf('Fig/Fig_%s/Binning', binStrategy); if isempty(dir(nameFigFolder_bin)), mkdir(nameFigFolder_bin), end
    
    load(sprintf('Data_compile/BinStrategy/BinInfo_n%d_nBins%d_%s.mat', nsubj, nBins, binStrategy))
    
    nSess_perBin_SD_allSubj = nan(nsubj, 1);
    dur_perBin_SD_allSubj = nSess_perBin_SD_allSubj;
    
    for isubj = 1:nsubj
        
        subjName = subjList{isubj};
        nblocks = nblocks_allSubj(isubj);
        nSess = nblocks/nLoc;
        nSess_perBin_suggested = round(nSess/nBins);
        
        date_perSess = date_perSess_allSubj{isubj};
        gap_diffFromLastSess = gap_diffFromLastSess_allSubj{isubj};
        nSess_perBin = nSess_perBin_allSubj{isubj};
        dur_perBin = dur_perBin_allSubj{isubj};
        gap_btwBin = gap_btwBin_allSubj{isubj};
        
        fprintf('\n============================')
        fprintf('\n%s nSess=%d, Total dur = %d days ', subjName, nSess, days_total_allSubj(isubj))
        fprintf('\n============================\n')
        
        fprintf('   Number of sessions per bin: [%s] SD=%.1f (suggest ~%d)\n', num2str(nSess_perBin'), std(nSess_perBin), nSess_perBin_suggested)
        fprintf('   Duration (in days) per bin: [%s] SD=%.1f\n', num2str(dur_perBin'), std(dur_perBin))
        fprintf('   Gaps (in days) between bins: [%s] (max=%d)\n', num2str(gap_btwBin'), max(gap_diffFromLastSess))
        
        nSess_perBin_SD_allSubj(isubj) = std(nSess_perBin);
        dur_perBin_SD_allSubj(isubj) = std(dur_perBin);
    end % isubj
    
    [nSess_perBin_SD_ave, ~, ~, nSess_perBin_SD_SEM] = getCI(nSess_perBin_SD_allSubj, 2, 1);
    [dur_perBin_SD_ave, ~, ~, dur_perBin_SD_SEM] = getCI(dur_perBin_SD_allSubj, 2, 1);
    
    fprintf('\n================ ACROSS SUBJ (n=%d) ===============\n', nsubj)
    fprintf('   SD of no. sessions per bin (ideally 0): ave=%.1f, SD=%.1f\n', nSess_perBin_SD_ave, nSess_perBin_SD_SEM)
    fprintf('   SD of duration (in days) per bin (ideally, 1): ave=%.1f, SD=%.1f', dur_perBin_SD_ave, dur_perBin_SD_SEM)
    fprintf('\n=================================================\n')
    
    % ====== plot ======
    
    % 1. number of sessions per bin
    nSess_perBin_allSubj_mtx = reshape(cell2mat(nSess_perBin_allSubj), nBins, nsubj)';
    basicFxn_drawBars(nSess_perBin_allSubj_mtx, [],colors_bin , {'Bin 1', 'Bin 2', 'Bin 3'}, [], [], 1, 0, sprintf('n=%d %s\nNumber of sessions per bin', nsubj, binStrategy), 0, [400 400]);
    saveas(gcf, sprintf('%s/n%d_nSessPerBin.jpg', nameFigFolder_bin, nsubj))
    
    % 2. duration (in days) per bin
    dur_perBin_allSubj_mtx = reshape(cell2mat(dur_perBin_allSubj), nBins, nsubj)';
    basicFxn_drawBars(dur_perBin_allSubj_mtx, [], colors_bin, {'Bin 1', 'Bin 2', 'Bin 3'}, [], [], 1, 0, sprintf('n=%d %s\nDuration (in days) per bin', nsubj, binStrategy), 0, [400 400]);
    saveas(gcf, sprintf('%s/n%d_durPerBin.jpg', nameFigFolder_bin, nsubj))
    
    % 3. gap (in days) between bins
    gap_btwBin_allSubj_mtx = reshape(cell2mat(gap_btwBin_allSubj), nBins-1, nsubj)';
    basicFxn_drawBars(gap_btwBin_allSubj_mtx, [], colors_bin, {'Bin 2-Bin 1', 'Bin 3-Bin 2'}, linspace(0, 300, 5), [], 1, 1, sprintf('n=%d %s\nGap (in days) between bin', nsubj, binStrategy), 0, [400 400]);
    saveas(gcf, sprintf('%s/n%d_gapBtwBin.jpg', nameFigFolder_bin, nsubj))
end % flag_binStrategy

%%
function [nSess_perBin, dur_perBin, gap_btwBin] = fxn_getBinInfo(nBins, nSess, date_perSess, iSess_start_allBins, gap_diffFromLastSess)

iSess_start_allBins_ = [1, iSess_start_allBins', nSess+1];
nSess_perBin = nan(nBins, 1);
dur_perBin = nSess_perBin;
for iBin=1:nBins
    date_end = date_perSess{iSess_start_allBins_(iBin+1)-1};
    date_start = date_perSess{iSess_start_allBins_(iBin)};
    dur_perBin(iBin) = ceil(days(date_end - date_start));
    nSess_perBin(iBin) = iSess_start_allBins_(iBin+1) - iSess_start_allBins_(iBin);
end
gap_btwBin = gap_diffFromLastSess(iSess_start_allBins-1);% +1 to match the index system of gap_diffFromLastSess
assert(sum(nSess_perBin) == nSess)
end
% % 1. extract the number of days each observer took, and the days between sessions
% % 2. decide how to bin data
% % 3. based on the binned data, calculate (1) the number of days each bin
% % took (2) the number of days between bins
% 
% clc, close all, warning off, format compact
% 
% addpath(genpath('Data'))
% addpath(genpath('fxn_analysis_RC_v2'))
% addpath(genpath('XueCarrasco_JN/code'))
% 
% %%
% clc,close all
% subjList = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA',  'DT', 'CS', 'DU', 'SR', 'RC'};
% nblocks_allSubj = [200, 240, 210, 240, 220, 220, 215, 220, 205, 210, 205, 205, 205, 195, 205];
% indSubj = 1:15;
% indSubj = [1:6, 8:14]; % no AS or RC
% indSubj = 1:11;
% subjList = subjList(indSubj);
% nblocks_allSubj = nblocks_allSubj(indSubj);
% 
% nsubj = length(subjList);
% 
% nLoc = 5; % each session contains 5 blocks, each testing one loc once; subj always finished multiples of 5 blocks per day
% nBins = 3;
% 
% % Extract date information
% % outputs:
% %       days_total_allSubj: nsubj x 1, total duration of data collection in days for each subj
% %       date_perSess_allSubj: nsubj x 1 / nSess x 1, the date (str, in datetime format) on which each session was completed
% %       days_diffFromLastSess_allSubj: nsubj x 1 / [nSess-1 x 1], the gap from last session in days
% 
% days_total_allSubj = nan(nsubj, 1);
% date_perSess_allSubj = cell(nsubj, 1);
% gap_diffFromLastSess_allSubj = cell(nsubj, 1);
% 
% for isubj = 1:nsubj
%     
%     subjName = subjList{isubj};
%     nblocks = nblocks_allSubj(isubj);
%     filesDir = dir(sprintf('Data/%s/%s_exp_B*', subjName, subjName));
%     if strcmp(subjName, 'SP'), filesDir = filesDir(81:end); end
%     assert(length(filesDir) == nblocks)
%     
%     nSess = nblocks/nLoc;
%     
%     % extract the date of each session
%     date_perSess = cell(nSess, 1);
%     gap_diffFromLastSess = nan(nSess-1, 1);
%     for iSess = 1:nSess
%         iblock = (iSess-1)*nLoc+1;
%         val_year = str2double(filesDir(iblock).name(15:18));
%         val_month = str2double(filesDir(iblock).name(19:20));
%         val_date= str2double(filesDir(iblock).name(21:22));
%         date_perSess{iSess} = datetime(val_year, val_month, val_date);
%         
%         if iSess>1
%             gap_diffFromLastSess(iSess-1) = round(days(date_perSess{iSess} -date_perSess{iSess-1}));
%         end
%     end % iSess
%     
%     date_perSess_allSubj{isubj} = date_perSess;
%     gap_diffFromLastSess_allSubj{isubj} = gap_diffFromLastSess;
%     
%     % Calculate the total duration (in days)
%     d1 = date_perSess{1};
%     dLast = date_perSess{end};
%     days_total_allSubj(isubj) = round(days(dLast-d1));
%     
% end
% 
% %% combine algorithm-based and manual search
% nSess_perBin_allSubj = cell(nsubj, 1);
% dur_perBin_allSubj{isubj} = nSess_perBin_allSubj;
% gap_btwBin_allSubj{isubj} = nSess_perBin_allSubj;
% iSess_start_perBin_allSubj = nSess_perBin_allSubj;
% iSess_end_perBin_allSubj = nSess_perBin_allSubj;
% 
% flag_binStrategy = input('       >>> 1=Algorithm-generated strategy; 2=Manually decided strategy, 0=equal: ');
% switch flag_binStrategy
%     case 0
%         binStrategy = 'equal';
%         nSess_perBin_allSubj = {[13, 13, 13]
%             [13, 13, 13]
%             [13, 13, 13]
%             [13, 13, 13]
%             [13, 13, 13]
%             [13, 13, 13]
%             [13, 13, 13]
%             [13, 13, 13]
%             [13, 13, 13]
%             [13, 13, 13]
%             [13, 13, 13]
%             [13, 13, 13]
%             [13, 13, 13]
%             [13, 13, 13]
%             [13, 13, 13]};
%         nSess_perBin_allSubj = nSess_perBin_allSubj(indSubj);
%         
%     case 1
%         binStrategy = 'algorithm';
%         nCombBest = 5;
%         nSess_perBin_min = 10; % the min. number of sessions each bin contains
%         
%     case 2
%         binStrategy = 'manual';
%         nSess_perBin_allSubj = {[12, 14, 14]
% [17, 15, 16]
% [14, 15, 13]
% [16, 17, 15]
% [13, 16, 15]
% [15, 17, 12]
% [14, 16, 13]
% [12, 16, 16]
% [15, 13, 13]
% [13, 15, 14]
% [13, 12, 16]}; %15-RC
%         nSess_perBin_allSubj = nSess_perBin_allSubj(indSubj);
% end
% 
% for isubj = 1:nsubj
%     clc
%     subjName = subjList{isubj};
%     nblocks = nblocks_allSubj(isubj);
%     nSess = nblocks/nLoc;
%     date_perSess = date_perSess_allSubj{isubj};
%     gap_diffFromLastSess = gap_diffFromLastSess_allSubj{isubj};
%     nSess_perBin_suggested = round(nSess/nBins);
%     
%     if flag_binStrategy==1
%         %%%%%% ALGORITHM-BASED %%%%%%
%         % the starting session of all but the first bin
%         iSess_start_allComb = combvec((nSess_perBin_min+1):nSess, (nSess_perBin_min+1)*2:nSess);
%         nComb = size(iSess_start_allComb, 2);
%         
%         val_allComb = nan(nComb, 1);
%         for iComb = 1:nComb
%             iSess_start_allBins = iSess_start_allComb(:, iComb);
%             if iSess_start_allBins(1)>=iSess_start_allBins(2)
%                 continue
%             else
%                 %------------------------%
%                 [nSess_perBin, dur_perBin, gap_btwBin] = fxn_getBinInfo(nBins, nSess, date_perSess, iSess_start_allBins, gap_diffFromLastSess);
%                 %------------------------%
%                 
%                 if any(nSess_perBin<nSess_perBin_min)
%                     continue
%                 else
%                     val_allComb(iComb) = std(nSess_perBin) - sum(gap_btwBin) + sum(dur_perBin);
%                 end
%             end
%         end % iComb
%         
%         [val_min, iComb_best] = sort(val_allComb);
%         
%         fprintf('\n============================')
%         fprintf('\n%s nSess=%d Total dur = %d days [max. gap=%d] (suggested nSess per bin): %d', ...
%             subjName, nSess, days_total_allSubj(isubj), max(gap_diffFromLastSess), nSess_perBin_suggested)
%         fprintf('\n============================\n')
%         
%         for iCombBest = 1:nCombBest
%             fprintf('\n*** [#%d] Val = %.2f***\n', iCombBest, val_min(iCombBest))
%             
%             %%%%%%%%%%%%
%             iSess_start_allBins = iSess_start_allComb(:, iComb_best(iCombBest));
%             [nSess_perBin, dur_perBin, gap_btwBin] = fxn_getBinInfo(nBins, nSess, date_perSess, iSess_start_allBins, gap_diffFromLastSess);
%             %%%%%%%%%%%%
%             
%             fprintf('Number of sessions per bin: [%s] SD=%.1f\n', num2str(nSess_perBin'), std(nSess_perBin))
%             fprintf('Gaps (in days) between bins: [%s]\n', num2str(gap_btwBin'))
%             fprintf('Duration (in days) per bin: [%s] \n', num2str(dur_perBin'))
%         end
%         iSess_start_allBins = iSess_start_allComb(:, iComb_best(1));
%         %-------------------------%
%         [nSess_perBin, dur_perBin, gap_btwBin] = fxn_getBinInfo(nBins, nSess, date_perSess, iSess_start_allBins, gap_diffFromLastSess);
%         %-------------------------%
%     else %%%%%% MANUAL %%%%%%
%         nSess_perBin = nSess_perBin_allSubj{isubj}.';
%         nSess = sum(nSess_perBin);
%         iSess_start_allBins = cumsum(nSess_perBin);
%         iSess_start_allBins = iSess_start_allBins(1:nBins-1)+1;
%         %-------------------------%
%         [nSess_perBin_, dur_perBin, gap_btwBin] = fxn_getBinInfo(nBins, nSess, date_perSess, iSess_start_allBins, gap_diffFromLastSess);
%         %-------------------------%
%         assert(sum(nSess_perBin_-nSess_perBin)==0)
%     end % if flagLoadBinStrategy
%     
%     nSess_perBin_allSubj{isubj} = nSess_perBin;
%     dur_perBin_allSubj{isubj} = dur_perBin;
%     gap_btwBin_allSubj{isubj} = gap_btwBin;
%     iSess_start_perBin_allSubj{isubj} = [1, cumsum(nSess_perBin(1:nBins-1)')+1];
%     iSess_end_perBin_allSubj{isubj} = cumsum(nSess_perBin');
% end % isubj
% 
% save(sprintf('Data_compile/BinStrategy/BinInfo_n%d_nBins%d_%s.mat', nsubj, nBins, binStrategy), '*_allSubj')
% 
% %% display the info of the best binning strategy (either by algorithm or by manual selection)
% clc
% load(sprintf('Data_compile/BinStrategy/BinInfo_n%d_nBins%d_%s.mat', nsubj, nBins, binStrategy))
% 
% nSess_perBin_SD_allSubj = nan(nsubj, 1);
% dur_perBin_SD_allSubj = nSess_perBin_SD_allSubj;
% 
% for isubj = 1:nsubj
%     
%     subjName = subjList{isubj};
%     nblocks = nblocks_allSubj(isubj);
%     nSess = nblocks/nLoc;
%     nSess_perBin_suggested = round(nSess/nBins);
%     
%     date_perSess = date_perSess_allSubj{isubj};
%     gap_diffFromLastSess = gap_diffFromLastSess_allSubj{isubj};
%     nSess_perBin = nSess_perBin_allSubj{isubj};
%     dur_perBin = dur_perBin_allSubj{isubj};
%     gap_btwBin = gap_btwBin_allSubj{isubj};
%     
%     fprintf('\n============================')
%     fprintf('\n%s nSess=%d, Total dur = %d days ', subjName, nSess, days_total_allSubj(isubj))
%     fprintf('\n============================\n')
%     
%     fprintf('   Number of sessions per bin: [%s] SD=%.1f (suggest ~%d)\n', num2str(nSess_perBin'), std(nSess_perBin), nSess_perBin_suggested)
%     fprintf('   Duration (in days) per bin: [%s] SD=%.1f\n', num2str(dur_perBin'), std(dur_perBin))
%     fprintf('   Gaps (in days) between bins: [%s] (max=%d)\n', num2str(gap_btwBin'), max(gap_diffFromLastSess))
%     
%     nSess_perBin_SD_allSubj(isubj) = std(nSess_perBin);
%     dur_perBin_SD_allSubj(isubj) = std(dur_perBin);
% end
% 
% [nSess_perBin_SD_ave, ~, ~, nSess_perBin_SD_SEM] = getCI(nSess_perBin_SD_allSubj, 2, 1);
% [dur_perBin_SD_ave, ~, ~, dur_perBin_SD_SEM] = getCI(dur_perBin_SD_allSubj, 2, 1);
% 
% fprintf('\n================ ACROSS SUBJ (n=%d) ===============\n', nsubj)
% fprintf('   SD of no. sessions per bin (ideally 0): ave=%.1f, SD=%.1f\n', nSess_perBin_SD_ave, nSess_perBin_SD_SEM)
% fprintf('   SD of duration (in days) per bin (ideally, 1): ave=%.1f, SD=%.1f', dur_perBin_SD_ave, dur_perBin_SD_SEM)
% fprintf('\n=================================================\n')
% 
% %% plot
% close all
% colors_bin = [0, 0, 0; .5, .5, .5; .75, .75, .75];
% 
% % 1. number of sessions per bin
% nSess_perBin_allSubj_mtx = reshape(cell2mat(nSess_perBin_allSubj), nBins, nsubj)';
% basicFxn_drawBars(nSess_perBin_allSubj_mtx, [],colors_bin , {'Bin 1', 'Bin 2', 'Bin 3'}, [], [], 1, 0, 'Number of sessions per bin', 0, [400 400]);
% 
% % 2. duration (in days) per bin
% dur_perBin_allSubj_mtx = reshape(cell2mat(dur_perBin_allSubj), nBins, nsubj)';
% basicFxn_drawBars(dur_perBin_allSubj_mtx, [], colors_bin, {'Bin 1', 'Bin 2', 'Bin 3'}, [], [], 1, 0, 'Duration (in days) per bin', 0, [400 400]);
% 
% % 3. gap (in days) between bins
% gap_btwBin_allSubj_mtx = reshape(cell2mat(gap_btwBin_allSubj), nBins-1, nsubj)';
% basicFxn_drawBars(gap_btwBin_allSubj_mtx, [], colors_bin, {'Bin 2-Bin 1', 'Bin 3-Bin 2'}, linspace(0, 300, 5), [], 1, 1, 'Gap (in days) between bin', 0, [400 400]);
% 
% %%
% function [nSess_perBin, dur_perBin, gap_btwBin] = fxn_getBinInfo(nBins, nSess, date_perSess, iSess_start_allBins, gap_diffFromLastSess)
% 
% iSess_start_allBins_ = [1, iSess_start_allBins', nSess+1];
% nSess_perBin = nan(nBins, 1);
% dur_perBin = nSess_perBin;
% for iBin=1:nBins
%     date_end = date_perSess{iSess_start_allBins_(iBin+1)-1};
%     date_start = date_perSess{iSess_start_allBins_(iBin)};
%     dur_perBin(iBin) = ceil(days(date_end - date_start));
%     nSess_perBin(iBin) = iSess_start_allBins_(iBin+1) - iSess_start_allBins_(iBin);
% end
% gap_btwBin = gap_diffFromLastSess(iSess_start_allBins-1);% +1 to match the index system of gap_diffFromLastSess
% assert(sum(nSess_perBin) == nSess)
% end