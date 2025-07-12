function noisePatch = exp_createNoisePatch(noiseParam, stimParam, ip)

% where is this function used???

filtered_noise = exp_CreateFilteredNoise(noiseParam);
noisePatch = noiseParam.ratio_base + filtered_noise * noiseParam.ratio_noiseInTgt;

noisePatch_pre3 = min(max(noisePatch .* stimParam.mask,0),1);
imshow(noisePatch_pre3), title(['patch #', num2str(ip)])
