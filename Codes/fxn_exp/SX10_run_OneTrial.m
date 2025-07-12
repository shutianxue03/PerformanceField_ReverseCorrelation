function run = SX10_run_OneTrial(exp_mode, EL, run, params, demoFlag)
% originally written by Antoine Barbot, adapted by Shutian Xue

screen = params.screen;
stim = params.stim;
time = params.time;
soundParams = params.soundParams;
window = screen.wPtr;
itrial = run.itrial;

%% present PH
exp_drawPH(1, stim, screen, run.iCuedLoc)
exp_drawFixCross(screen,stim)
time1 = Screen('Flip', window);

if run.ON, rd_eyeLink('trialstart', window, {EL, run.itrial, screen.centerX, screen.centerY, screen.rad});
else, WaitSecs(time.fix);
end
fixation = 1;

%% check fixation
stopThisTrial = 0;

while (GetSecs < time1 + time.fix) && (~stopThisTrial)
    WaitSecs(.01);
    if run.ON, fixation = rd_eyeLink('fixcheck', window, {screen.centerX, screen.centerY, screen.rad});end
    if ~fixation
        DrawFormattedText(params.screen.wPtr, sprintf('Please fixate'),'center', 'center');
        Screen('Flip', screen.wPtr); WaitSecs(1);
        stopThisTrial = 1;
    end
end
run.fixation = fixation;
if ~run.fixation, return,end

if (run.iblock * itrial == 1) && (exp_mode == 2) && demoFlag, KbWait;end

%% present PH
exp_drawPH(1, stim, screen,run.iCuedLoc)
exp_drawFixCross(screen,stim)
Screen('Flip', window);
WaitSecs(.3); % to ensure the presence of a stable image when go back to exp, after the eyetracking is broken
if run.iCuedLoc > 1, exp_drawFixCross(screen,stim), end

%% present PH & stim
% fprintf('tgtPrs = %d.\n', run.tgtPrs)
exp_prsStim(stim, screen, run.stimPatch_all) % stim were created in SX9
exp_drawPH(1, stim, screen, run.iCuedLoc)
time2 = Screen('Flip', window);  % display the stimuli

while (GetSecs < time2 + time.stim) && (~stopThisTrial)
    WaitSecs(.01);
    if run.ON, fixation = rd_eyeLink('fixcheck', window, {screen.centerX, screen.centerY, screen.rad});end
    if ~fixation
        DrawFormattedText(params.screen.wPtr, sprintf('Please fixate'),'center', 'center');
        Screen('Flip', screen.wPtr); WaitSecs(1);
        stopThisTrial = 1;
    end
end
run.fixation = fixation;
if ~run.fixation;return,end

if ~run.ON, fixation = 1;end % if not using eyetracker, just reset fixation
run.fixation = fixation;

if (run.iblock * itrial == 1) && (exp_mode == 2) && demoFlag, KbWait;end

%% wait for response
% exp_drawFixCross(screen, stim)
exp_drawPH(1,stim,screen, run.iCuedLoc); Screen('Flip', screen.wPtr);
WaitSecs(time.wait);

makeBeep(soundParams, soundParams.freqNeutral, time.feedback)

%% detect key
switch fixation
    case 0 % not fixating
        if run.ON,Eyelink('command','draw_text 100 100 42 Fixation break');end
    case 1 % check keypress
        run.RT= eps;
        while run.RT<.001
            run = exp_getKey(soundParams, time, stim, run);
            if run.RT < .001 % the key is pressed before the beep ends
                run.check = 0;
                DrawFormattedText(params.screen.wPtr, 'Please press the key after the beep.','center', 'center'); 
                Screen('Flip', screen.wPtr); 
                WaitSecs(.5);
                exp_drawPH(1, stim, screen, run.iCuedLoc)
                Screen('Flip', window);
                return
            end
        end
        run = exp_checkResp(params, run,1); 
        if run.ON, Eyelink('message', 'EVENT_ClearScreen');end
end

if (run.iblock * itrial == 1) && (exp_mode == 2) && demoFlag, KbWait; end

%% ITI (no fixation cross)
exp_drawPH(1, stim, screen, run.iCuedLoc)
Screen('Flip', window);
WaitSecs(time.ITI);

if (run.iblock * itrial == 1) && (exp_mode == 2) && demoFlag, KbWait; end

