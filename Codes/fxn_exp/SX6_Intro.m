function SX6_Intro(params)

limit0to1 = @(x) min(max(x,0),1);

stim = params.stim;
noise = params.noise;
screen  = params.screen;
soundParams = params.soundParams;
time = params.time;

%% write the line
DrawFormattedText(screen.wPtr, stim.line_instr ,'center', 'center', screen.black);

%% draw the sample stimuli
gabor = exp_CreateGabor(stim, .8, 0);
filtered_noise = exp_CreateFilteredNoise(noise);

noisyGabor = limit0to1(noise.ratio_base + gabor*noise.ratio_gaborInTgt + filtered_noise).*stim.mask + noise.ratio_base*(1-stim.mask);
noisePure = limit0to1(noise.ratio_base + filtered_noise).*stim.mask + noise.ratio_base*(1 - stim.mask);

% display 2 types of stimuli
for n = 1:2
    if n==1, image =  noisyGabor; else, image = noisePure;  end
    textureLoc = Screen('MakeTexture',screen.wPtr,image,[],[],1);
    Screen('DrawTexture', screen.wPtr, textureLoc, [], CenterRectOnPointd(stim.rect, screen.centerX + (-1)^n*stim.ecc*1, screen.centerY+30))
end

%% display the gabor
gabor_ = limit0to1(noise.ratio_base + gabor * noise.ratio_gaborInTgt).*stim.mask + noise.ratio_base*(1 - stim.mask);
textureLoc = Screen('MakeTexture',screen.wPtr, gabor_,[],[],1);
Screen('DrawTexture', screen.wPtr, textureLoc, [], CenterRectOnPointd(stim.rect, screen.centerX - stim.ecc*2, screen.centerY+30))
DrawFormattedText(screen.wPtr, '(          )' , screen.centerX - 442, screen.centerY + 15, screen.black);

% show them
Screen('Flip', params.screen.wPtr); 
KbStrokeWait;

