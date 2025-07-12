function run = EL_preExp(EL,params,run)
design = params.design;
screen = params.screen;

%%
Eyelink('command','clear_screen');
if run.itrial == 1 % do calibration
    calib_result = EyelinkDoTrackerSetup(EL);
    if calib_result == EL.TERMINATE_KEY, return, end
end

Eyelink('message', 'STARTING TRIALS');

if  (Eyelink('isconnected') == EL.notconnected), return, end % cancel if eyeLink is not connected
Eyelink('command', 'record_status_message ''Block %d of %d, Trial %d of %d''', run.iblock, design.nBlocks, run.itrial, design.nTrialsPerBlock);
Eyelink('message', 'Trial ID: %d', run.itrial);

ncheck = 0; % number of fixation checks
status.fixBreak    = 0;  % status of fixation (0: not fixating at the center)
status.recording = 0; % status of recording

while ~status.fixBreak % exit the loop once start to fixate or start to record
    if ~status.recording
        Eyelink('startrecording');	% start recording
        WaitSecs(.1); % start recording 100 msec before just to be safe
        
        key=1; while key, key = EyelinkGetKey(EL);  end  % dump any pending local keyss
        
        status.recording = Eyelink('checkrecording'); 	% Check if we are recording: if not, report an error.
        if ~status.recording
            status.recording = 1;
            Eyelink('message', 'RECORD_START');
        else
            status.recording = 0;
            Eyelink('message', 'RECORD_FAILURE');
        end
        
    else % already recording
        Eyelink('command','clear_screen 0');
        status.fixBreak  = exp_checkBreak(screen, EL);
        ncheck = ncheck + 1;
        calib_result   = EyelinkDoTrackerSetup(EL);
        if calib_result == EL.TERMINATE_KEY, return, end
        status.recording = 0;
    end
end
run.ON = 1;
run.fixBreak = status.fixBreak;
run.ncheck = ncheck;
Eyelink('message', 'TRIAL_START %d', run.itrial);
Eyelink('message', 'SYNCTIME');		% zero-plot time for EDFVIEW