
%% save data of last run
% save(participant.nameFile_lastRun ,'real_sequence','sequence','stimulus','response','timing','params','visual','scr','design');

%% display thank-you message
rubber([]);
Screen('FillRect', scr.main, visual.bgColor, scr.rect); Screen('Flip', scr.main); WaitSecs(1);
DrawFormattedText(scr.main,'Thanks, you have finished this part of the experiment.','center', 'center', visual.black);
Screen('Flip',scr.main);

%% turn off scr
WaitSecs(2);
Screen('CloseAll');

%% turn off EL recording (after each block)
if constant.EYETRACK
    if ~exist(el.eyeDataDir,'dir'), mkdir(el.eyeDataDir); end
    rd_eyeLink('eyestop', scr.main, {el.eyeFile, el.eyeDataDir});
end

%% rename eyedata file
if constant.EYETRACK
    string_datetime = real_sequence.string_datetime;
    string_b = real_sequence.string_b;
    eyeDataFileName = sprintf('eyedata/%s_E%d_%s_%s.edf', ...
        participant.subjName, constant.expMode, string_b, string_datetime);
    movefile('eyedata/xx.edf', eyeDataFileName)
end

%% save params of this subj
paramsMode = {'stair', 'prac', 'conStim'};
if find(constant.expMode == [1,3])
    rmfield(params, 'UD');
    save(sprintf('%s/Params_%s.mat', constant.nameFolder, paramsMode{constant.expMode}), ...
        'params', 'constant' , 'design', 'participant', 'scr', 'stimulus', 'timing', 'visual', 'scr')
    fprintf('Saved\n')
else
    fprintf('\n\nPractice mode: Params not saved.\n')
end

%% turn off sound track
for ff = [800, 600, 400], makeBeep(ff, .2), end
PsychPortAudio('Close', params.pahandle);

%%
fprintf(1,'\n\nDuration: %.1f min\n\n.', (toc)/60);

