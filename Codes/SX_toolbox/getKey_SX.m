function [keyPress, time] = getKey_SX(keyPress)

global keys constant

% use 'findKeyPressed' to test/find the key being pressed

[KeyIsDown, t, keyCode] = KbCheck;
time = [];

if KeyIsDown
    key = find(keyCode==1);
    if sum(key == keys.respButtons(:))
%         keyPress = find(key == keys.respButtons); % key assignment is set in 'getKeyAssignment'
        keyPress = key;
        time    = t;

    elseif key == keys.stopKey
        if constant.EYETRACK, Eyelink('Shutdown');end
        makeBeep(1e3, .2)
        makeBeep(1e3, .2)
        Screen('CloseAll'); error('ALERT: Experiment stopped by the user!');
        
    end
end