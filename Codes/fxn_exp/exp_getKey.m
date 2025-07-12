function run = exp_getKey(soundParams, time, stim, run)
% Description
% This script detect and record the key pressed
% NOTE: if unsure about the key code of a key, use 'checkKeyCode' to check it.
% KbName('UnifyKeyNames')

% KbName('UnifyKeyNames');

CorrecKeyPressed = 0;
t_start = GetSecs;

while ~CorrecKeyPressed
    [KeyIsDown, t_end, keyCode] = KbCheck;
    if KeyIsDown, keyPressed = find(keyCode);
        if sum(keyPressed == stim.keysCodes)
            CorrecKeyPressed = 1;
            run.key = keyPressed;
            run.RT = t_end - t_start;
            run.check = CorrecKeyPressed;
        elseif keyPressed == stim.stopKey
            %             if keyCode(KbName('q'))
            if run.ON, Eyelink('Shutdown');end
            makeBeep(soundParams, 1000, time.feedback)
            makeBeep(soundParams, 1000, time.feedback)
            Screen('CloseAll'); error('ALERT: Experiment stopped by the user!'); 
        end
    end
end
