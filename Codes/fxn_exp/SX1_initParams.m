
function params = SX1_initParams(screen_mode,exp_mode)

if exp_mode < 4 % not simulation
    screen =  SX2_initScreen(screen_mode);
    ppd = screen.ppd; % just to simplify
else, ppd = 32; screen.ratio_base = .5;
end

%% exp design
design.nPrsLoc = 5; % no. of stimuli presented
design.targetOrNot = [0,1]; % detection task, present or absent
ratio_prs = .5; % prop. of tgt-prs trials
design.nrep = 2; % number of times a stim is repeated

if exp_mode < 4 % not simulation
    design.ratio_prs = ratio_prs;
    
    if exp_mode == 2, loc_list = [1,3,2,4,5]; design.nCuedLoc = design.nPrsLoc;
    else, loc_list = randperm(design.nPrsLoc); 
        design.nCuedLoc = length(loc_list); 
    end
    design.loc_list = loc_list;
    
    switch exp_mode
        case 0, disp('Titration') 
            design.nTrialsPerStairPerLoc = 50; % default 50, must be multiples of 2
            design.nStairs = 2;
        case 1, disp('Experiment') 
            design.nTrialsPerStairPerLoc = 100; % default 100, must be multiples of 4
            design.nStairs = 1;
        case 2, disp('Practice')
            design.nTrialsPerStairPerLoc = 10; % default 20
            design.nStairs = 1;
        case 3, disp('Quick titration')
            design.nTrialsPerStairPerLoc = 10; % default 50
            design.nStairs = 1;
    end
    design.nBlocks = length(loc_list);
    
    if exp_mode == 2 % practice (data not saved)
        design.ind = allcomb(1:design.nStairs, design.targetOrNot);
        design.trialInd = repmat(design.ind, design.nTrialsPerStairPerLoc * 5,1);
        design.nTrials = size(design.trialInd,1);
        design.nTrialsPerBlock = design.nTrials/design.nBlocks;
        locInd = repmat(loc_list, 2*design.nTrialsPerStairPerLoc, 1);
        design.trialInd = [locInd(:), design.trialInd(randperm(design.nTrials), :)];
    else % not practice, data saved
        ntrials = design.nTrialsPerStairPerLoc * design.nStairs * length(loc_list);% all blocks
        design.nTrials = ntrials;
        trialInd_stair_tgt = []; % which loc is cued | which stair | target or not
        for nL = 1: length(loc_list)
            rep = ntrials/length(loc_list)/design.nStairs;
            idx_stair1 = ones(rep,1);
            idx_targetOrNot_stair1 = [ones((round(design.nTrialsPerStairPerLoc) * ratio_prs), 1); zeros(round(design.nTrialsPerStairPerLoc * (1-ratio_prs)), 1)];
            idx_stair2 = ones(rep,1)*2;
            idx_targetOrNot_stair2 = idx_targetOrNot_stair1;
            
            idx_both_1 = [idx_stair1, idx_targetOrNot_stair1];
            idx_both_1_rand = idx_both_1(randperm(rep),:);
            idx_both_2 = [idx_stair2, idx_targetOrNot_stair2];
            idx_both_2_rand = idx_both_2(randperm(rep),:);
            
            % experiment, need only ONE staircases
            if exp_mode, idx_both_2_rand = []; end
            trialInd_stair_tgt = [trialInd_stair_tgt; [idx_both_1_rand; idx_both_2_rand]];
        end
        
        locInd = repmat(loc_list, ntrials/length(loc_list), 1);
        design.trialInd = [locInd(:), trialInd_stair_tgt];
        if exp_mode == 1, design.trialInd = []; end
        design.nTrialsPerBlock = ntrials/design.nBlocks;
    end
    
    clc
    fprintf('==============\n%d blocks x %d trials =% d trials.\n', design.nBlocks, design.nTrialsPerBlock, design.nTrials)
    disp('Location tested:')
    loc_list
    disp('==============')
    
    %% time (ms)
    time.fix = .1; % fixation time, in sec, default = .2
    time.ITI = .1;  % only PH, default = .2
    time.stim = .1; % .stimulus presentation time, in sec
    time.wait = .3; % default = .5, can start to respond after this amount of time after the stim offset,
    time.rest = 5; % mandatory rest time between blocks
    time.feedback = .1; % duration of the feedback beep, in sec
    
end % below will be run in simulation

%% stim and sizes
stim.ppd = ppd; % produced by SX2
if exp_mode < 4, stim.noncue_color = screen.black; end
stim.stim5 = 1; % 1 = present all 5 stimuli; 0 = only present 1
stim.PH_loc_box = allcomb([1,-1],[1,-1]);
stim.ecc = 6 * ppd; % eccentricity
stim.ecc_box = [0,0;-stim.ecc,0; 0,-stim.ecc;stim.ecc,0;0,+stim.ecc]; % each row 1: center, LHM, UVM, (RHM, LVM
stim.ecc_box = stim.ecc_box(1:design.nPrsLoc,:);
% if size(stim.ecc_box,1) ~= design.nPrsLoc, makeBeep(.1, [1000 1000]), error('ALERT: Number of locations unmatched.'),end
stim.PH_dist = 2 * ppd;
stim.PH_leng = .75 * ppd;
stim.PH_sz = stim.PH_dist / 2 + stim.PH_leng; % in dva
if screen_mode == 0, stim.PH_wid = .04 * ppd; else, stim.PH_wid = .1 * ppd;  end
stim.fixSz = .3 * ppd; % length of fixation arms
gaborSF =  2;
stim.gaborSF = gaborSF;

stim.gaborSD = .8; % AF
stim.gabor_sz = 3; % the width of stim, or object sz
stim.aper_sz = 3; % The visible part of the stimulus, or image size
stim.sin_sz = .8; % AF
stim.sinpower = 5; % AF's default

stim.gabor_psz = round(stim.gabor_sz * ppd); 
stim.aper_psz = round(stim.aper_sz * ppd); 
stim.sin_psz = round(stim.sin_sz * ppd); % AF

stim.mask  = exp_CreateCircularApertureSin(stim);
% stim.mask = ones(stim.aper_psz, stim.aper_psz);
stim.rect = [0, 0, stim.aper_psz, stim.aper_psz];

stim.text_sz = round(5*sqrt(ppd)); % in dva
stim.cue_color = 180; % not that bright
stim.targetOri = 90; % in deg, horizontal

%% noise
noise.noiseCST = 0.2;
noise.ratio_gaborInTgt = .5; % in Ian's code: .5
noise.ratio_noiseInTgt = 1;
noise.ratio_base = screen.ratio_base; % set by the bgColor in SX2_initScreen
noise.SF_multiplier = 2;
noise.sz = stim.aper_sz;
noise.psz = round(noise.sz * ppd);
noise.SF_low = gaborSF/noise.SF_multiplier;
noise.SF_high = gaborSF*noise.SF_multiplier;
noise.filterSF_all_log = linspace(log2(noise.SF_low), log2(noise.SF_high), 15);
noise.filterSF_all = 2.^(noise.filterSF_all_log); % on a log scale
noise.filterSF_all_ppd = noise.filterSF_all/ppd;

noise.ori_base = 0;
noise.ori_diff = 180;
noise.fix_contrast = 1; %1 = the rmsContrast of the noise is fixed across trials; 0 = allowing more energy fluctuation across trials
noise.ori_low = [];  % noise.ori_base-(noise.ori_diff/2);
noise.ori_high = []; % noise.ori_base+(noise.ori_diff/2);
noise.ppd = ppd;
noise.scr = screen;

% rotateNoise % the function that rotates a pre-made noise patch; not in use

%% staircase
if exp_mode < 4 % not simulation
    if exp_mode == 2, nAlpha = 5; stairParams.initThresh = [3, .6]; % easier for practice to help observer to learn the task
    else, nAlpha = 500; stairParams.initThresh = [.6, .1]; 
    end
    stairParams.PF = @arbWeibull; % weibull fxn with arbitrary performance level at which threshold is help
    % created by Michael,
    % used to be @PAL_CumulativeNormal;
    stairParams.whichStair = 2; % 1 = best Pest (mode); 2=PEST (mean)
    stairParams.threshPerformance = .70; %   the perf at which the performance is held
    stairParams.lastPosterior = [];
    stairParams.alphaRange = 10.^(linspace(log10(stairParams.initThresh(2)), log10(stairParams.initThresh(1)), nAlpha));
    stairParams.fitBeta = 2;
    stairParams.fitGamma = 0.5;
    stairParams.fitLambda = 0.01;
    
    if exp_mode == 0 % premapping the alpha to explore a wide range fo alpha space
        stairParams.updateAfterTrial = 10;
        stairParams.preUpdateLevels = linspace(stairParams.initThresh(1), stairParams.initThresh(2), stairParams.updateAfterTrial);
    else
        stairParams.updateAfterTrial = 1;
        stairParams.preUpdateLevels = 0;
    end
    % stairParams.stopRule = design.nTrialsPerStairPerLoc;
    % stairParams.startStepsize = .01;
    % stairParams.minStepsize = 10^-4;
    
    %% keys
    systemName = computer;
    if systemName(1) == 'M' % Mac
        stim.keysCodes = [9,13]; % F and J; windows: 89 98; mac: 9,13; Linux: 42,45
        stim.stopKey = 20; % for q; mac: 20, Linux:25
        stim.space = 44; % space key, linux = 66, Mac=44
    elseif systemName(1) == 'G' % linux, at least for the L1 computer
        stim.keysCodes = [42,45]; % F and J; windows: 89 98; mac: 9,13; Linux: 42,45
        stim.stopKey = 25; % for q; mac: 20, Linux:25
        stim.space = 66; % space key, linux = 66, Mac=44
    else
        disp('system not identified.')
    end
    
    %% lines
    exp_modes = {'Titration','Experiment','Practice', 'Quick titration', 'simulation'};
    stim.line_instr2 = '\n\n\nWelcome! Your task is to report whether the gabor is present or not.\n\n\nF = YES    J = NO\n\n\nPress any key to start.\n';
    stim.line_instr1 = exp_modes{exp_mode+1};
    stim.line_instr = strcat(stim.line_instr1,stim.line_instr2);
    
    stim.line_end = 'You have finished all blocks.\n\nThank you for your participation!';
    Screen('TextSize', screen.wPtr, stim.text_sz);
    
    %% sound
    InitializePsychSound(1); % 1 for precise timing
    % Open the  audio device
    pahandle = PsychPortAudio('Open');
    device = PsychPortAudio('GetDevices');
    soundParams.fs = device.DefaultSampleRate;
    soundParams.pahandle =  pahandle;
    soundParams.freqCorrect = 1000;
    soundParams.freqNeutral = 800;
    soundParams.freqWrong = 600;
    soundParams.amp = 1;
    
    %% organize
    params.screen = screen;
    params.design = design;
    params.stairParams = stairParams;
    params.time = time;
    params.soundParams = soundParams;
end

params.noise = noise;
params.stim = stim;

if exp_mode < 4, for ff = [600, 800, 1000], makeBeep(soundParams, ff, time.feedback*2), end, end

