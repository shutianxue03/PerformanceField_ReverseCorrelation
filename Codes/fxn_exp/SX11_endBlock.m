function record = SX11_endBlock(exp_mode, fileName, params, run, record)
screen = params.screen;
design = params.design;
iblock = run.iblock;
iblockAll = record.iblock_current-1;

%% display performance
% iblock_current = mod(run.iblock,5); if iblock_current == 0, iblock_current = 5;end
accPerBlock = nanmean(record.correctness(iblockAll,1:design.nTrialsPerBlock)) ;
answ_ = record.answer(iblockAll, 1:design.nTrialsPerBlock);
prs_ = record.tgtPrs(iblockAll, 1:design.nTrialsPerBlock);
nHit = sum(answ_ & prs_);
nFA = sum(answ_ & (1-prs_));
[dprime, criterion] =  SX_sim06_SDT(nHit, nFA, design.nTrialsPerBlock/2);
    
fprintf('B%d/%d L%d: accuracy = %.2f, dprime = %.2f, criterion = %.2f, duration: %d secs.\n', iblock, design.nBlocks, design.loc_list(iblock), accPerBlock, dprime, criterion, round(record.duration(iblock)));

% DrawFormattedText(screen.wPtr, 'Saving data.\n' ,'center', 'center', screen.black); Screen('Flip', screen.wPtr);
record.accPerBlock(iblock) = accPerBlock;
record.dprime(iblock) = dprime;
record.criterion(iblock) = criterion;

%% display accuracy
% line_report = sprintf('You have finished block %d/%d.\nAccuracy = %d%%\n',run.iblock,design.nBlocks,round(accPerBlock *100));
line_report = sprintf('accuracy = %d%%\n', round(accPerBlock *100));
if run.iblock ~= design.nBlocks % not the last block
    for irest = 1:params.time.rest
        DrawFormattedText(screen.wPtr, sprintf('%s\nNext block will be ready in %d secs.\n',line_report,params.time.rest-irest+1) ,'center', 'center', screen.black);
        Screen('Flip', screen.wPtr);
        WaitSecs(1);
    end
else % the last block
    DrawFormattedText(screen.wPtr, line_report,'center', 'center', screen.black);
    Screen('Flip', screen.wPtr);
    WaitSecs(1);
end

clear accPerBlock ans line_instru run

%% save data by block (except practice)
fileName.time = datestr(now,'yyyymmddTHHMM');
if exp_mode ~= 2
    if record.iblock_current <= 10 % 1-9 % record.iblock_current is always 1 bgger than the real iblock!!
        fileName.formal = sprintf('%s_B00%dL%d_%s',fileName.brief, record.iblock_current-1, record.iCuedLoc(record.iblock_current-1,1), fileName.time);
    elseif (record.iblock_current >= 11) && (record.iblock_current<=100) % 10-99
        fileName.formal = sprintf('%s_B0%dL%d_%s',fileName.brief, record.iblock_current-1, record.iCuedLoc(record.iblock_current-1,1), fileName.time);
    else % 100 and above
        fileName.formal = sprintf('%s_B%dL%d_%s',fileName.brief, record.iblock_current-1, record.iCuedLoc(record.iblock_current-1,1), fileName.time);
    end
    if exp_mode ~= 1, record.run_allTrials = []; end  % if not exp blocks, do not save stimuli
    
    save(fileName.formal, 'record', 'params', 'exp_mode', 'fileName')
end
record.fileName = fileName;


