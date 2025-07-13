
%% save data by block (except practice)
if constant.expMode~=2
    if constant.expMode==1
        iblock_toSave = run.iblock; 
    else
        iblock_toSave = params.lastBlock;
    end
    string_datetime = datestr(now,'yyyymmdd_HHMM');
    if iblock_toSave<10, string_b = sprintf('b0%d', iblock_toSave);
    else, string_b = sprintf('b%d', iblock_toSave);
    end
    nameFile_perBlk = sprintf('%s/%s_E%d_%s_%s', constant.nameFolder, participant.subjName, constant.expMode, string_b, string_datetime);
%     dirFile_perBlk = sprintf('%s/%s', participant.nameFolder, nameFile_perBlk); %strcat('Data/',datName,'/',nameFile_perBlk);
    save(nameFile_perBlk, 'participant','real_sequence','sequence','stimulus','response','timing','params','visual','scr','design');
    real_sequence.string_datetime = string_datetime;
    real_sequence.string_b = string_b;
end

%% display accuracy on the screen and countdown
acc_perBlock = mean(real_sequence.iscor);
line_report = sprintf('Accuracy = %d%%\n', round(acc_perBlock*100));
time_rest = 5;
if iblock_local ~= design.nBlocksPerSession % not the last block
    for irest = 1:time_rest
        DrawFormattedText(scr.main, sprintf('%s\nNext block will be ready in %d secs.\n',line_report, time_rest-irest+1) ,'center', 'center', visual.black); Screen('Flip', scr.main);
        WaitSecs(1);
    end
else % the last block
    DrawFormattedText(scr.main, line_report,'center', 'center', visual.black); Screen('Flip', scr.main);
    WaitSecs(1);
end

%% display performance in the command window
if constant.expMode ~= 2
    oriResp = real_sequence.oriResp;
    stimOri = real_sequence.stimOri;
    iscor = real_sequence.iscor;
    pHit = mean(iscor(stimOri==1));
    pFA = mean(1-iscor(stimOri==2));
    [dprime, criterion] =  SX_SDT(pHit, pFA);
    
    fprintf('\n\nB%d/%d (B%d/%d): accuracy = %d%%, dprime = %.2f, criterion = %.2f (neg. prefer LEFT), confidence=%.1f\n\n\n', ...
        iblock_local, design.nBlocksPerSession, iblock, design.nBlockTotal, ...
        round(acc_perBlock*100), dprime, criterion, mean(real_sequence.confidence));
    
    % DrawFormattedText(screen.wPtr, 'Saving data.\n' ,'center', 'center', screen.black); Screen('Flip', screen.wPtr);
    real_sequence.acc_perBlock = acc_perBlock;
    real_sequence.dprime = dprime;
    real_sequence.criterion = criterion;
end
