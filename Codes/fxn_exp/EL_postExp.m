function  [continue_flag ,run] = EL_postExp(EL, run,params)

if run.ON, Eyelink('message', 'TRIAL_END %d', run.itrial);Eyelink('Stop recording');end

if run.fixBreak == 1
    continue_flag = 0;
    DrawFormattedText(params.screen.wPtr, 'Please fixate.','center', 'center', params.screen.black); Screen('Flip', params.screen.wPtr);
    if ~EL.ON; fprintf(1,' Fixation break - trial added, now total of %i trials',params.design.nTrialsPerBlock); end
elseif run.fixBreak==0 
    continue_flag = 1;
end