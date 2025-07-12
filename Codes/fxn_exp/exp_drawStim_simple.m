
function [stim_tgt, stim_noise] = exp_drawStim_simple(params, cst)

stim = params.stim;
noise = params.noise;

noise.ratio_base = .5;

limit0to1 = @(x) min(max(x,0),1);

% generate noisy gabor and noise
[gabor, ~] = exp_CreateGabor(stim, cst, 0);
filtered_noise = exp_CreateFilteredNoise(noise);
stim_tgt = limit0to1(noise.ratio_base + gabor*noise.ratio_gaborInTgt + filtered_noise).*stim.mask + noise.ratio_base*(1-stim.mask);
stim_noise = limit0to1(noise.ratio_base + filtered_noise).*stim.mask + noise.ratio_base*(1 - stim.mask);





