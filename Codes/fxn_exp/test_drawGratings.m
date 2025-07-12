function run = test_drawGratings(screen,run,design,stim,noise)
stim.phase = rand*2*pi; % phase at each trial is randomized
run.contrast = run.threshold;
[run.gabor, run.gaborCarrier, run.gaborModulator] = exp_CreateRaisedGrating(screen, stim.gabor_sz, stim.gratinSF, stim.flatSpread, stim.phase, 0, stim.targetOri, plotQ);

for nGrating = 1:design.numLoc
    %  the texture of Grating
    % maniputlate: noise.SF_low, noise.SF_high, noise.ori_low, noise.ori_high
    run.filtered_noise = exp_CreateFilteredNoise(screen, noise.size, noise.SF_low, noise.SF_high, noise.ori_low, noise.ori_high, noise.contrast(nGrating), noise.fix_contrast);
    run.noisyGrating = .5 + run.gabor*.5  + run.filtered_noise;
    
    run.Grating = min(max(run.Grating,0),1);
    run.texture = cat(3,run.Grating,stim.mask);
    
    % the loc of Grating
    switch nGrating
        case 1, run.GratingRect = [];% 1=center ([] is the default loc)
        case 2, run.GratingRect = CenterRectOnPointd(stim.rect,screen.centerPix(1) - stim.ecc, screen.centerPix(2)); % 2= LH
        case 3, run.GratingRect = CenterRectOnPointd(stim.rect,screen.centerPix(1), screen.centerPix(2)-stim.ecc); % 3=UVM
    end
    run.GratingTexture = Screen('MakeTexture',screen.wPtr,run.texture,[],[],1);
    Screen('DrawTexture', screen.wPtr, run.GratingTexture, [], run.GratingRect)
end

exp_drawPH(1,stim,screen,run), Screen('Flip', screen.wPtr);
WaitSecs(stim.t_stim);

% placeholders + response cue
exp_drawPH(1,stim,screen,run), Screen('Flip', screen.wPtr);