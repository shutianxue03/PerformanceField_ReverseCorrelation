function [] = A54_FFTP(expt,X,sm)
% Generates X noise images and plot of them containing the
% energy in Fourier space for all spatial frequencies and orientations in
% the noise
% expt = screen and stimulus parameter
% X = number of simulated patches per subject
% sm = smooth 1 or 0
%% generate new noise
scr = expt.scr;
stiPar=expt.stiPar;
nSubs = 6;
nTrials = X;
nSims = (nTrials)*nSubs;
%mask  = A05_CreateCircularApertureSin(scr, stiPar.gaborsiz, stiPar.sinsiz, [], stiPar.apersiz);

for ii = 1:nSubs
    %ns = mkNoise(scr, stiPar.apersiz, stiPar.noise_lf, stiPar.noise_hf, stiPar.noise_contrast, 1,sm).*mask;
    ns = mkNoise(scr, stiPar.apersiz, stiPar.noise_lf, stiPar.noise_hf, stiPar.noise_contrast, 1,sm);
    %     ns = 0.5 + ns;
    %     ns = min(max(ns,0),1);
    patch{ii} = ns;
    clear ns
end
%% Fourier transform of the noise across all stim
for ii = 1:length(patch)
    F = fft2(patch{ii});
    xF = mean(F,1);
    x(:,ii)=abs(xF);
end
%% plot
F_x = mean(x,2); F_x(1)=0;
g_size = 2; % grating size
bar([0:(length(xF)-1)]/g_size,F_x,'FaceAlpha',0.7,'EdgeColor',[1 1 1]); af_darkTheme;
xlabel('cycles/deg');xlim([0 8]);


%% function to make the noise
    function [filterednoise] = mkNoise(scr, noisesiz, noise_lf,noise_hf, contrast, fixcontrast,sm)
        
        noisepsiz = A15_angle2pix(scr, noisesiz);
        sampling_period = 1/A15B_pix2angle(scr,1);
        nyquist         = sampling_period/2;
        
        freq     = linspace(0,nyquist,floor(noisepsiz/2)+1);
        faxis    = [ freq  fliplr(freq(2:ceil(noisepsiz/2))) ];
        faxis    = fftshift(faxis);
        [fgrid_x,fgrid_y] = meshgrid(faxis,faxis);
        [ogrid, fgrid]    = cart2pol(fgrid_x,fgrid_y);
        
        %% Draw filter for frequency domain
        fFilter = double(fgrid > noise_lf & fgrid < noise_hf);
        
        if sm
            smooth = normpdf(1:nyquist,mean(1:nyquist),nyquist/10); % was set to 10 in experiment
            fFilter = conv2(fFilter,smooth,'same');
        end
        
        %% Generate Gaussian noise and apply the filters
        noise = randn(noisepsiz, noisepsiz);
        fn    = fftshift(fft2(noise));
        filterednoise = real(ifft2(ifftshift(1.*fFilter.*fn)));
        
        if fixcontrast ==0
            
            filterednoise = filterednoise*contrast - mean(filterednoise(:));
        elseif fixcontrast == 1
            
            filterednoise = filterednoise/std(filterednoise(:))*contrast - mean(filterednoise(:));
        end
        
    end
end


