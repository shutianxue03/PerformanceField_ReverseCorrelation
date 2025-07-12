function [filterednoise, fFilter] = exp_CreateFilteredNoise(noise)

nyquist  = noise.ppd/2;
freq     = linspace(0,nyquist,floor(noise.psz/2)+1);
faxis    = [freq, fliplr(freq(2:ceil(noise.psz/2))) ];
faxis    = fftshift(faxis);
[fgrid_x,fgrid_y] = meshgrid(faxis,faxis);
[ogrid, fgrid]    = cart2pol(fgrid_x,fgrid_y); 

%% Draw filter for frequency domain
fFilter = double(fgrid > noise.SF_low & fgrid < noise.SF_high);
% smooth = normpdf(1:nyquist,mean(1:nyquist), nyquist/10); % nyquist/10 is used by AF
% fFilter = conv2(fFilter,smooth,'same'); %added this

%% Generate Gaussian noise and apply the filters
noise_patch = randn(noise.psz, noise.psz);
fn    = fftshift(fft2(noise_patch));
filterednoise = real(ifft2(ifftshift(1.*fFilter.*fn))); %changed this

if noise.fix_contrast == 0 % % rmsContrast of the noise is NOT fixed across trials to allow more energy fluctuation
    filterednoise = filterednoise*noise.noiseCST - mean(filterednoise(:));
else % rmsContrast of the noise is fixed across trials
    filterednoise = filterednoise/std(filterednoise(:))*noise.noiseCST - mean(filterednoise(:));
end


