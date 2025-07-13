function startBlock_SX(iblock_local, design, el)

global params visual scr keys constant
if constant.expMode ~= 2
    DrawFormattedText(scr.main, sprintf('%s\nBlock %d/%d\n',...
        constant.names_expModes{constant.expMode}, iblock_local, design.nBlocksPerSession), 'center', 'center', visual.black);
    DrawFormattedText(scr.main, sprintf('Press space to calibrate\nPress Q to quit\nPress any other keys to continue\n'), scr.centerX  - 200, scr.centerY + 50, visual.black);
else % practice
    DrawFormattedText(scr.main, sprintf('%s\nBlock %d/%d\n',...
        constant.names_expModes{constant.expMode}, iblock_local, design.nBlocksPerSession), 'center', 'center', visual.black);
    DrawFormattedText(scr.main, sprintf('Press Q to quit\nPress any other keys to continue\n'), scr.centerX  - 200, scr.centerY + 50, visual.black);
end
Screen('Flip', scr.main);

%% wait for key presses
CorrecKeyPressed = 0; % unmute
% CorrecKeyPressed = 1;keyIsDown=1;keyPressed=44; % mute
while ~ CorrecKeyPressed
    [keyIsDown, ~, keyCode] = KbCheck; % -1: gather all input devices); % the input should NOT be the screen ind!!
    if keyIsDown && (sum(keyCode) > 0), keyPressed = find(keyCode); CorrecKeyPressed = 1; end
end

switch keyPressed
    case keys.space % Press space to calibrate
        if constant.EYETRACK
            [~, exitFlag] = rd_eyeLink('calibrate', scr.main, el);
            if exitFlag, return, end
        end
    case keys.stopKey % Press Q to quit
        if constant.EYETRACK, Eyelink('Shutdown');end
        makeBeep(1e3, .2), makeBeep(1e3, .2)
        Screen('CloseAll'); sca; error('ALERT: Experiment stopped by the user!');
end

Screen(scr.main,'FillRect',visual.bgColor,[]);
Screen(scr.main, 'Flip');
