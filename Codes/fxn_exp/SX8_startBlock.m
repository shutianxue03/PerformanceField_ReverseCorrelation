function SX8_startBlock(params, EL, iblock_current)

locNames = {'CENTER', 'LEFT', 'UPPER', 'RIGHT', 'LOWER'};

scr = params.screen; 
stim = params.stim; 
design = params.design;
loc_list = design.loc_list; % the list of loc ind to be run in this session
nLoc_tested = length(loc_list);

% determine the index of the current block within the 5-blocks bundle
iblock_current_mod = mod(iblock_current, 5);
if iblock_current_mod == 0, iblock_current_mod = nLoc_tested;end
ib_current_name = locNames{loc_list(iblock_current_mod)};
DrawFormattedText(scr.wPtr, sprintf('Block %d\n\n The tested location is %s\n', iblock_current, ib_current_name), 'center', 'center', scr.black);
DrawFormattedText(scr.wPtr, sprintf('Press space to calibrate\nPress Q to quit\nPress any other keys to continue\n'), scr.centerX  - 200, scr.centerY + 100, scr.black);

Screen('Flip', scr.wPtr);

%% wait for key presses
CorrecKeyPressed = 0;

while ~ CorrecKeyPressed
    [~, ~, keyCode] = KbCheck;
    if sum(keyCode) > 0, keyPressed = find(keyCode); CorrecKeyPressed = 1; end
end

switch keyPressed
    case stim.space % Press space to calibrate
        if EL.ON
            [~, exitFlag] = rd_eyeLink('calibrate', scr.wPtr, EL);
            if exitFlag, return, end
        end
        
    case stim.stopKey % Press Q to quit
        if EL.ON, Eyelink('Shutdown');end
        makeBeep(params.soundParams, 1000, params.time.feedback)
        Screen('CloseAll'); sca; error('ALERT: Experiment stopped by the user!');
end

