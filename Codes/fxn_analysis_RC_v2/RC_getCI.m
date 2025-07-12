

% compute classification image based on the noise field ONLY
% CI is generated for tgt-prs, abs and both trials, respectively

noisePrs = noise_both_perComb_p{iLoc}(1:ntrials*length(iLoc_all));
noiseAbs = noise_both_perComb_p{iLoc}(ntrials*length(iLoc_all)+1:end);

idxHit = boolean(data_both_allT_perComb{iLoc, 2}(1:ntrials*length(iLoc_all))); % column 2 is the resp (1=yes, 0=no)
idxFA = boolean(data_both_allT_perComb{iLoc, 2}(ntrials*length(iLoc_all)+1:end));
idxMiss = boolean(1-idxHit);
idxCR = boolean(1-idxFA);

% hit (prs and say yes)
noiseHit = cell2Mat3(noisePrs(idxHit));

% FA (abs and say yes)
if sum(idxFA) == ntrials, noiseFA = [];
else, noiseFA = cell2Mat3(noiseAbs(idxFA));
end

% miss (prs and say no)
if sum(idxMiss) == ntrials, noiseMiss = [];
else, noiseMiss = cell2Mat3(noisePrs(idxMiss));
end

% CR (absent and say no)
noiseCR = cell2Mat3(noiseAbs(idxCR));

fprintf('%s: pHit=%.2f, pFA=%.2f, pMiss=%.2f, pCR=%.2f\n', namesLocComb{iLoc}, mean(idxHit), mean(idxFA), mean(idxMiss), mean(idxCR))

%% get CI
CI_mean_ = mean(noiseHit,3) + mean(noiseFA,3) - mean(noiseCR,3) - mean(noiseMiss,3);
CI_var_ = var(noiseHit,[],3) + var(noiseFA,[],3) - var(noiseCR,[],3) - var(noiseMiss,[],3);

CI_mean{iLoc} = CI_mean_;
CI_var{iLoc} = CI_var_;

clear noisePrs noiseAbs noiseHit noiseFA noiseMiss noiseCR



