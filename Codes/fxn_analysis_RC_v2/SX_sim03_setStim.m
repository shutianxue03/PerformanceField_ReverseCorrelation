

limit0to1 = @(x) min(max(x,0),1);
% rmsCST = @(patch, sz)sqrt(sum(sum((patch-mean(patch(:))).^2))/(sz^2)); %the common meaning of 'cst'

gaborPhase = nan(1,ntrials);

patchPrs = cell(1,ntrials);
patchAbs = patchPrs;
noisePrs = patchPrs;
noiseAbs = patchPrs;
noisePrs_mask = patchPrs;
noiseAbs_mask = patchPrs;

for itrial = 1:ntrials
    [gabor, phase] = exp_CreateGabor(stim, stim.gaborCST,0);
    gaborPhase(itrial) = phase;
    
    filtered_noise1 = exp_CreateFilteredNoise(noise);
    filtered_noise2 = exp_CreateFilteredNoise(noise);
    
    noisePrs{itrial} = filtered_noise1;
    noiseAbs{itrial} = filtered_noise2;
    
    noisePrs_mask{itrial} = limit0to1(noise.ratio_base+filtered_noise1).*stim.mask + noise.ratio_base*(1-stim.mask);
    noiseAbs_mask{itrial} = limit0to1(noise.ratio_base+filtered_noise2).*stim.mask + noise.ratio_base*(1-stim.mask);
    
    patchPrs{itrial} = limit0to1(noise.ratio_base+gabor*noise.ratio_gaborInTgt+filtered_noise1).*stim.mask + noise.ratio_base*(1-stim.mask);
    patchAbs{itrial} = limit0to1(noise.ratio_base+filtered_noise2).*stim.mask + noise.ratio_base*(1-stim.mask);
end

noise_both = cat(2, noisePrs, noiseAbs);
noise_mask_both = cat(2, noisePrs_mask, noiseAbs_mask);
stim_both = cat(2,patchPrs, patchAbs);



