    

 
% rd_eyeTrackTester.m
%
% This is designed to be a very simple "experiment" that can be used to
% test and/or demonstrate the use of rd_eyeLink.m.

subjectID = 'test';
eyeDataDir = 'eyedata';

eyeFile = sprintf('%s  1111%s', subjectID, datestr(now, 'mmdd'));

nTrials = 5;
rad = 70; % radius of allowable eye movement in pixels
% makeBeep(.1, [600,800,1000])

%% Screen
screenNumber = max(Screen('Screens'));
% [window rect] = Screen('OpenWindow', screenNumber);
[window rect] = Screen('OpenWindow', screenNumber, [], [0 0 800 600]); % for testing
[cx cy] = RectCenter(rect);
Screen('TextSize', window, 24);
Screen('TextColor', window, 255);
Screen('TextFont', window, 'Verdana');
Screen('FillRect', window, 128)
Screen('Flip', window);

%% Initialize eye tracker
[el exitFlag] = rd_eyeLink('eyestart', window, eyeFile);
if exitFlag, return, end

%% Calibrate eye tracker
[cal exitFlag] = rd_eyeLink('calibrate', window, el);
if exitFlag, return, end

%% Present trials
for iTrial = 1:nTrials
    
    fprintf('\n\nTrial %d', iTrial)
    
    %% instruction
    Screen('FillRect', window, 128)
    DrawFormattedText(window, 'Press any key to start', 'center', 'center', 0);
    Screen('Flip', window);
    KbWait;
    
    %% present fixation cross
    DrawFormattedText(window, '+', 'center', 'center', 255);
    timeFix = Screen('Flip', window);
    
    %% check fixation
    % when the eyetracker is recording and the subject is holding fixation
    rd_eyeLink('trialstart', window, {el, iTrial, cx, cy, rad});
    
    %% present first stimulus  
    DrawFormattedText(window, '+', 'center', 'center', 255);
    DrawFormattedText(window, 'STIM 1', cx-200, 'center');
    timeStim1(iTrial) = Screen('Flip', window);
    Eyelink('Message', 'EVENT_STIM1');
    
    stopThisTrial = 0;
    while GetSecs < timeStim1(iTrial) + 0.45 && ~stopThisTrial
        WaitSecs(.01);
        fixation = rd_eyeLink('fixcheck', window, {cx, cy, rad});
        if ~fixation
            fprintf('\nBROKE FIXATION! (stim 1)')
            Beeper('high')
            stopThisTrial = 1;
        end
        % alt: fixationBreakTasks(fixation)
    end
    fix1(iTrial) = fixation;
    
    if stopThisTrial, continue, end
    
    %% present second stimulus
    DrawFormattedText(window, '+', 'center', 'center');
    DrawFormattedText(window, 'STIM 2', cx+200, 'center');
    timeStim2(iTrial) = Screen('Flip', window, timeStim1(iTrial) + 0.5);
    Eyelink('Message', 'EVENT_STIM2');
    
    while GetSecs < timeStim2(iTrial) + 0.45 && ~stopThisTrial
        WaitSecs(.01);
        fixation = rd_eyeLink('fixcheck', window, {cx, cy, rad});
        if ~fixation
            fprintf('\nBROKE FIXATION! (stim 2)')
            Beeper('high')
            stopThisTrial = 1;
        end
    end
    fix2(iTrial) = fixation;
    
    Eyelink('Message', 'TRIAL_END');
end

%% Final blank
Screen('FillRect', window, 128)
Screen('Flip', window);

 %% Save the eye data and shut down the eye tracker
if ~exist(eyeDataDir,'dir'), mkdir(eyeDataDir), end
rd_eyeLink('eyestop', window, {eyeFile, eyeDataDir});

%% Close screen
makeBeep(.1, [1000,800,600])
Screen('CloseAll')

