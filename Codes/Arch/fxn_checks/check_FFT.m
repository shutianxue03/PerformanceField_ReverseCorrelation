
%% 
% ==============================
% test the patches used in my exp
% ==============================

clc
limit0to1 = @(x) min(max(x,0),1);
load('Data/SX/SXparams1')
scr = params.screen;
stim = params.stim;
noise = params.noise;
psz = stim.aper_psz;
sz = stim.aper_sz;
% ratio_base = .5;

prs = 0;
cstGabor = 1;
ntrials = 1000;%2*1e3;%*1e3;

dimInd = 2; % average 2D power to 1D power

% create mask
mask  = exp_CreateCircularApertureSin(stim);
for ratio_base = [.3, .4, .5] %.3:.1:.6 % higher, brighter

    noisePure_all = cell(ntrials,1);
    for itrial = 1:ntrials
        noisePure_all{itrial} = exp_CreateFilteredNoise(noise);
    end
    
    FFT2D_gabor = nan(ntrials, psz, psz);
    FFT2D_noisePure = FFT2D_gabor;
    FFT2D_tgt_nonMask = FFT2D_gabor;
    FFT2D_tgt_mask = FFT2D_gabor;
    
    FFT1D_gabor = nan(ntrials, psz);
    FFT1D_noisePure = FFT1D_gabor;
    FFT1D_tgt_nonMask = FFT1D_gabor;
    FFT1D_tgt_mask = FFT1D_gabor;
    
    for itrial = 1:ntrials
        % create Gabor
        [gabor, ~] = exp_CreateGabor(stim, cstGabor, 0); % random phase
        [FFT2D, FFT1D] = getFFT(gabor, dimInd);
        FFT2D_gabor(itrial, :, :) = FFT2D;
        FFT1D_gabor(itrial, :) = FFT1D;
        
        % Create pure noise
        %         noisePure = exp_CreateFilteredNoise(noise);
        noisePure = noisePure_all{itrial};
        [FFT2D, FFT1D] = getFFT(noisePure, dimInd);
        FFT2D_noisePure(itrial, :, :) = FFT2D;
        FFT1D_noisePure(itrial, :) = FFT1D;
        
        % targte without the circular mask
        if prs == 1, tgt = ratio_base + gabor*ratio_base + noisePure;
        else , tgt = ratio_base + noisePure;
        end
        tgt_nonMask = limit0to1(tgt);
        [FFT2D, FFT1D] = getFFT(tgt_nonMask, dimInd);
        FFT2D_tgt_nonMask(itrial, :, :) = FFT2D;
        FFT1D_tgt_nonMask(itrial, :) = FFT1D;
        
        % mask added
        tgt_mask = tgt_nonMask.*mask + .5 * (1-mask);
        [FFT2D, FFT1D] = getFFT(tgt_mask, dimInd);
        FFT2D_tgt_mask(itrial, :, :) = FFT2D;
        FFT1D_tgt_mask(itrial, :, :) = FFT1D;
        
    end
    
    % ================
    % take average
    FFT2D_gabor_ave = squeeze(mean(FFT2D_gabor, 1));
    FFT2D_noisePure_ave = squeeze(mean(FFT2D_noisePure, 1));
    FFT2D_tgt_nonMask_ave = squeeze(mean(FFT2D_tgt_nonMask, 1));
    FFT2D_tgt_mask_ave =squeeze(mean(FFT2D_tgt_mask, 1));
    
    FFT1D_gabor_ave = mean(FFT1D_gabor);
    FFT1D_noisePure_ave = mean(FFT1D_noisePure);
    FFT1D_tgt_nonMask_ave = mean(FFT1D_tgt_nonMask);
    FFT1D_tgt_mask_ave = mean(FFT1D_tgt_mask);
    % ================
    % plot
    figure('position', [0 0 600 800])
    subplot(4,3,1), imshow(gabor), axis square, title('Gabor')
    plotFFT(FFT2D_gabor_ave, FFT1D_gabor_ave, sz, 2)
    
    subplot(4,3,4), imshow(noisePure), axis square, title('Pure noise')
    plotFFT(FFT2D_noisePure_ave, FFT1D_noisePure_ave, sz, 5)
    
    subplot(4,3,7), imshow(tgt_nonMask), axis square, title('Target unmasked')
    plotFFT(FFT2D_tgt_nonMask_ave, FFT1D_tgt_nonMask_ave, sz, 8)
    
    subplot(4,3,10), imshow(tgt_mask), axis square, title('Target masked')
    plotFFT(FFT2D_tgt_mask_ave, FFT1D_tgt_mask_ave, sz, 11)
    
    set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
    sgtitle(sprintf('[smoothed] prs = %d, pedestal = %.1f', prs, ratio_base), 'FontSize',25)
end

%%
% figure('Position', [0 0 500 300])
% subplot(1,2,1), stem((1:sz*10)/sz, FFT1D_tgt_nonMask_ave(2:(sz*10+1))), axis square, xticks([1,2,4,5]), xline(2, 'r'); xline(1, 'g'); xline(4, 'g'); xlabel('SF (cpd)'), title('Unsmoothed, Unmasked')
% subplot(1,2,2), stem((1:sz*10)/sz, FFT1D_tgt_mask_ave(2:(sz*10+1))),       axis square, xticks([1,2,4,5]), xline(2, 'r'); xline(1, 'g'); xline(4, 'g'); xlabel('SF (cpd)'), title('Unmoothed, Masked')
% set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
% sgtitle(sprintf('Aperture size = %.1f deg', sz), 'FontSize',20)


%% 
% ==============================
% check AF's patches used for demo
% ==============================

% addpath(genpath('AF_test/'))
% load AF_test/testData_full.mat
% load exptParams.mat
% patch_AF = patch; clear patch
% stiPar = expt.stiPar;
% scr = expt.scr;
% 
% nAllTrials = size(trialsMat, 1);
% cueLocation = trialsMat(:,5);
% tgtType = nan(nAllTrials,1);
% for ii = 1:nAllTrials
%     if cueLocation(ii) == 1 || cueLocation(ii) == 3, tgtType(ii,:) = trialsMat(ii,3);
%     elseif cueLocation(ii) == 2 || cueLocation(ii) == 4, tgtType(ii,:) = trialsMat(ii,4);
%     end
% end
% cst_all = trialsMat(:,9); % contrast, 2 levels
% cst = unique(cst_all);
% 
% patch_prs_cst1 = patch_AF(tgtType == 1 & cst_all == cst(1));
% patch_prs_cst2 = patch_AF(tgtType == 1 & cst_all == cst(2));
% patch_abs_cst1 = patch_AF(tgtType == 0 & cst_all == cst(1));
% patch_abs_cst2 = patch_AF(tgtType == 0 & cst_all == cst(2));
% 
% %%
% for n=4 % plot the 2D and 1D power of 4 types of trials
%     switch n
%         case 1, pp =  patch_prs_cst1; sgtitle_ = 'prs cst1';
%         case 2, pp =  patch_prs_cst2;sgtitle_ = 'prs cst2';
%         case 3, pp =  patch_abs_cst1;sgtitle_ = 'abs cst1';
%         case 4, pp =  patch_abs_cst2;sgtitle_ = 'abs cst2';
%     end
%     ntrials = length(pp);
%     pp_fft2_2D = nan(ntrials, stiPar.gaborpsiz, stiPar.gaborpsiz);
%     pp_fft2_1D = nan(ntrials, stiPar.gaborpsiz);
%     
%     margDim = 2; % 1=average vertically, 2= average horizontally
%     % given the power plot of prs trials (clear spike at SF=2), margDim=1 seems to be correct
%     
%     for it = 1:ntrials
%         [fft2D, fft1D] = getFFT_(pp{it}, margDim, 1);
%         pp_fft2_2D(it, :, :) = fft2D;
%         pp_fft2_1D(it, :) = fft1D;
%     end
%     
%     plotFFT_(pp_fft2_2D, margDim)
%     sgtitle(sgtitle_, 'fontsize', 20)
%     clear pp
% end

%% helper fxns
function [fft2D, fft1D] = getFFT(pp, dimInd)
fft2D= abs(fft2(pp));
fft1D = abs(mean(fft2(pp), dimInd));
% fft1D = mean(abs(fft2(pp)),dimInd);
end

function plotFFT(ave_2D, ave_1D, sz, iplot)

subplot(4,3,iplot)
imagesc(fftshift(ave_2D(2:end, 2:end))), axis square
colorbar, colormap gray

subplot(4,3,iplot+1)
stem(log2((1:(sz*10))/sz), ave_1D(2:(sz*10+1))), axis square
xticks([0,1,2]),xticklabels([1,2,4])
xline(1, 'r');
xline(0, 'g');
xline(2, 'g');
xlabel('SF (cpd)')
ylabel('power')
end

