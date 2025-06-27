
% If you use this method or found this tutorial helpful please cite:
% Fernández, A., Okun, S., & Carrasco, M. (2022). Differential effects of endogenous and exogenous attention on sensory tuning. Journal of Neuroscience, 42(7), 1316-1327.
% Xue, S., Fernández, A., & Carrasco, M. (2024). Featural Representation and Internal Noise Underlie the Eccentricity Effect in Contrast Sensitivity. The Journal of Neuroscience: The Official Journal of the Society for Neuroscience, 44(3). doi: 10.1523/JNEUROSCI.0743-23.2023

noise.ppd = 32;
noise.psz = 3*noise.ppd ;
noise.SF_low = 1/noise.ppd;
noise.SF_high = 4/noise.ppd;
noise.fix_contrast = 1;
noise.noiseCST = 20/100;

stim.mask = 
[filtered_noise, fFilter] = exp_CreateFilteredNoise(noise);

limit0to1 = @(x) min(max(x,0),1);

stim_noise = limit0to1(.5 + filtered_noise).*stim.mask + .5*(1 - stim.mask);


figure
imagesc(filterednoise)
axis square

function [filterednoise, fFilter] = exp_CreateFilteredNoise(noise)

nyquist  = noise.ppd/2;
freq     = linspace(0,nyquist,floor(noise.psz/2)+1);
faxis    = [freq, fliplr(freq(2:ceil(noise.psz/2))) ];
faxis    = fftshift(faxis);
[fgrid_x,fgrid_y] = meshgrid(faxis,faxis);
[~, fgrid]    = cart2pol(fgrid_x,fgrid_y);

%% Draw filter for frequency domain
fFilter = double(fgrid > noise.SF_low & fgrid < noise.SF_high);
smooth = normpdf(1:nyquist,mean(1:nyquist), nyquist/10); % nyquist/10 is used by AF
fFilter = conv2(fFilter,smooth,'same'); %added this

%% Generate Gaussian noise and apply the filters
noise_patch = randn(noise.psz, noise.psz);
fn    = fftshift(fft2(noise_patch));
filterednoise = real(ifft2(ifftshift(1.*fFilter.*fn))); %changed this

if noise.fix_contrast == 0 % % rmsContrast of the noise is NOT fixed across trials to allow more energy fluctuation
    filterednoise = filterednoise*noise.noiseCST - mean(filterednoise(:));
else % rmsContrast of the noise is fixed across trials
    filterednoise = filterednoise/std(filterednoise(:))*noise.noiseCST - mean(filterednoise(:));
end

end
