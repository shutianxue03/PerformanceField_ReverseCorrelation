function exp_prsStim(stim, screen, stimPatch_all)

for iLoc = 1:length(stimPatch_all)
    % decide the loc of stimulus
    switch iLoc
        case 1, gaborRect = CenterRectOnPointd(stim.rect,screen.centerX, screen.centerY); % 2= LH;% 1=center ([] is the default loc)
        case 2, gaborRect = CenterRectOnPointd(stim.rect,screen.centerX - stim.ecc, screen.centerY); % 2= LH
        case 3, gaborRect = CenterRectOnPointd(stim.rect,screen.centerX, screen.centerY-stim.ecc); % 3=UVM
        case 4, gaborRect = CenterRectOnPointd(stim.rect,screen.centerX+ stim.ecc, screen.centerY); % 4=right
        case 5, gaborRect = CenterRectOnPointd(stim.rect,screen.centerX, screen.centerY+stim.ecc); % 5=LVM
    end
    textureLoc = Screen('MakeTexture', screen.wPtr, stimPatch_all{iLoc}, [], [], 1);
    Screen('DrawTexture', screen.wPtr, textureLoc, [], gaborRect)
end