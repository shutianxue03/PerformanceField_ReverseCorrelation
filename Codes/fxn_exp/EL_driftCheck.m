function eye = EL_driftCheck(eye,screen,run,expStruct)

driftCorr = 1;
if eye.mode && (bb > 1)
    driftCorr = EyelinkDoDriftCorrect(eye.el, screen.centerPix(1), screen.centerPix(2), 1, 1);
end
if eye.mode && ~driftCorr, EyelinkDoDriftCorrect(eye.el, 'c'); end

blank = NaN(1,expStruct.nTrials);
eye.results{run.nBlock} = struct('RT',blank,'key', blank,'correct',blank);

% breakFix{bb} = struct('check', {[]}, 'track', {[]}, 'count', 0, 'recent', 0); %REPLICATE IN OTHER EXP SCRIPTS
BF.count = 0; 
BF.recent = 0; 
BF.recalCount = 0;
eye.breakFix{run.nBlock} = BF;
