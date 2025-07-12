% simulate patches to test the correlation between energy bins
% Using AF's code and params
% could create white noise by switching CreateFilteredNoise_temp to CreateFilteredNoise_noBound

%% define params
phase = 1; % arbitrary
ori = 90; % arbitrary
cst_signal = .12;
cst_noise = .18;

stiPar_AF = expt.stiPar;
scr_AF = expt.scr;
sz = stiPar_AF.apersiz;
g_std = stiPar_AF.gaborstd; % 0.8
sf = stiPar_AF.gaborf; % 2

%% step 1: Gabor
gabor_AF = AF_mkGabor(scr_AF, sz, g_std, sf, phase, cst_signal, ori);
% here, the 2nd input 'sz' was originally 'gaborsiz' in AF's code, which takes the value 1.9
% However, the size of the output is 64, and only when sz = apersiz (2.1) the patch size is 70

%% step 2: masks
mask_AF  = AF_mkMask(scr_AF, sz, stiPar_AF.sinsiz, [], sz);

%% step 3: create multiple PRS/ABS patches
fft_gabor_AF_2D = abs(fft2(gabor_AF)); 
ntrials = 1.5e3;
nfft = round(sz*5);
patch0_AF_allT = cell(ntrials, 1);
patch1_AF_allT = patch0_AF_allT;
fft_noise_AF_1D = nan(ntrials, nfft);

limit0to1 = @(x) min(max(x,0),1);

plotFlag = 0;

if plotFlag, figure('Position', [2000 0 2000 2000]), end
for itrial = 1:ntrials
    
    fre_high = stiPar_AF.noise_hf;
    fre_low = stiPar_AF.noise_lf;
        [noise_AF, fFilter_AF] = CreateFilteredNoise_temp(scr_AF, sz, fre_low, fre_high, cst_noise, 1);
%     noise_AF = CreateFilteredNoise_noBound(scr_AF, sz, cst_noise, 1);
    
    patch0_AF_allT{itrial} = limit0to1(.5 + noise_AF).*mask_AF + .5*(1 - mask_AF);
    patch1_AF_allT{itrial} = limit0to1(.5 + gabor_AF*.5 + noise_AF).*mask_AF + .5*(1-mask_AF);
    fft_noise_AF_2D = abs(fft2(noise_AF)); 
    fft_noise_AF_1D(itrial, :) = fft_noise_AF_2D(2:nfft+1, 2);
    
    % figure
    if plotFlag
        subplot(4,4,3), imshow(gabor_AF)
        subplot(4,4,4), imshow(fft_gabor_AF_2D)
        subplot(4,4,7), imshow(noise_AF)
        subplot(4,4,8), imshow(fft_noise_AF_2D)
        
        subplot(4,4,15), imshow(patch1_AF_allT{itrial}), title('[PRS] AF')
        subplot(4,4,16), imshow(patch0_AF_allT{itrial}), title('[ABS] AF')
        
        waitforbuttonpress
    end
end
fprintf('DONE\n')

%% get energy
clc
patch_all = patch0_AF_allT;
lumiBG = patch_all{1}(1,1);
ORI_all = -90:7.5:90; nORI = length(ORI_all);
nSF = 25; SF_log_all = linspace(log2(1), log2(4), nSF); SF_all = 2.^SF_log_all;

e3D_sim = nan(ntrials, nORI, nSF);
for iORI = 1:nORI
    for iSF = 1:nSF
        SF = SF_all(iSF);
        ORI = ORI_all(iORI);
        filter1 = CreateGabor(scr_AF, sz, g_std, SF, 0, 1, ORI);
        filter2 = CreateGabor(scr_AF, sz, g_std, SF, pi/2, 1, ORI);
%         filter1 = filter1.*mask_AF;
%         filter2 = filter2.*mask_AF;
        filter1 = filter1/(sum(filter1(:).^2));
        filter2 = filter2/(sum(filter2(:).^2));
        
        for itrial = 1:ntrials
            tempImg = patch_all{itrial}-lumiBG;
            sinImg = tempImg(:)'*filter1(:);
            cosImg = tempImg(:)'*filter2(:);
            e3D_sim(itrial, iORI, iSF) = sqrt(sinImg^2 + cosImg^2);
        end
    end
end


fprintf('energy DONE\n')
% xaxis_tuning = {ORI_all, SF_all};
xaxis_tuning = {ORI_all, SF_log_all};
xticks_contour_tuning = {[1, 7, 13, 19, 25], [1, 7, 13, 19, 25]};
xticks_tuning = {-90:45:90; 0:.5;2};


%% helper fxn - noise without bound
function filterednoise = CreateFilteredNoise_noBound(scr, noisesiz, contrast, fixcontrast)

noisepsiz = angle2pix(scr, noisesiz);
noise = randn(noisepsiz, noisepsiz);

if fixcontrast ==0
    filterednoise = noise*contrast - mean(noise(:));
elseif fixcontrast == 1
    filterednoise = noise/std(noise(:))*contrast - mean(noise(:));
end
end
