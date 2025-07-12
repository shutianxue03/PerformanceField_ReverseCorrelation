
% This script simiulates an ideal observer that detect the presence of a
% Gabor (2 cpd, horizontally orientated, contrast undefined) embedded in
% bandpass noise (contains SF from 1-4 cpd, contrast undefined)

% The tempate is a prewhitened Gabor given that the pixels in the noise
% patch are correlated.

% On each trial, the response variable (RV) is calculated by taking the dot
% product of the fourier component of the prewhitened template
% (template_prw) and the prewhitened stimuli: prewhitened noise (noise_prw)
% or the prewhitened noisy gabor (gabor_prw+noise_prw). Then, the RV is
% compared to the optimal criterion, which is half the distance between the
% distribution of noise (N) and noisy Gabor (S+N).

% Theoretically, an ideal observer makes the decision by comparing the
% stimulus distribution to the criterion. The center of the noise
% distribution (mu_stimN) is 0; The center of the noisy Gabor distribution
% (mu_stimG) is the dot product of the prewhitened Gabor and prewhitened
% template (the same Gabor). Two distributions share the same width.

clc
close all

%% settings
mode_prw = 1; % 0=using white noise; 1=using bandpass noise and prewhiten stimuli and template;
nTrials = 1e2; % number of simulated trials
% define Gabor contrast
nGaborCST = 30; gaborCST_log_all = linspace(-3, 0, nGaborCST); gaborCST_ln_all = 10.^gaborCST_log_all;
% define noise contrast
nNoiseCST = 5; noiseCST_all = linspace(.2, 1, nNoiseCST);
cst_ln_template = 100/100; % the contrast of the template
gaborSF = 2; % Gabor SF
SF_lb = 1; SF_ub = 4; % lower and upper bound of the bandpass filter
gaborSD = .8; % Gaussian

ppd = 32; % pixel per degree
aper_sz = 3; % stim diameter, in dva
aper_psz = aper_sz*ppd; % stim diameter, in pixel
nyquist  = ppd/2;
targetOri = 90; % 90 = horizontal

stim.aper_sz=aper_sz;
stim.aper_psz=aper_psz;
stim.targetOri=targetOri;
stim.phase=0;
stim.gaborSF=gaborSF;
stim.gaborSD = gaborSD;
stim.gabor_sz = aper_sz;

perf_thresh = .7;

%% create the bandpass filter in the Fourier domain
freq = linspace(0, nyquist, floor(aper_psz/2)+1);
faxis = [freq, fliplr(freq(2:ceil(aper_psz/2))) ];
[fgrid_x, fgrid_y] = meshgrid(faxis, faxis);
[~, fgrid]    = cart2pol(fgrid_x, fgrid_y);
filter_fft = double(fgrid >= SF_lb & fgrid <= SF_ub);

%% create & prewhiten the template
phase_template = 0;
template = exp_CreateGabor(stim, cst_ln_template);
template_prw = fxn_prewhiten(template, filter_fft);

%% calculate the predicted center of N
mu_theo_N = 0;

%%
sigma_sim_allN = nan(1, nNoiseCST);
sigma_theo_allN = sigma_sim_allN;
thresh_ln_allN = sigma_sim_allN;

clc

for iNoiseCST = 1:nNoiseCST
    
    noiseCST = noiseCST_all(iNoiseCST);
    fprintf('Noise#%d/%d (%.3f)...\n', iNoiseCST, nNoiseCST, noiseCST)
    
    % calculate the sigma for both gabor-present (S+N) and absent distributions (N)
    if mode_prw % if bandpass noise is used
        %=========================================%
%         sigma_theo = sqrt(sum(abs(template_prw).^2))*noiseCST;
        sigma_theo = sqrt(sum(real(template_prw).^2 + imag(template_prw).^2))*noiseCST;
        sigma_theo = sqrt(sum(real(template_prw).^2) + sum(imag(template_prw).^2))*noiseCST;
        %=========================================%
    else
        sigma_theo = sqrt(sum(template.^2, 'all'))*noiseCST;
    end
    sigma_theo_allN(iNoiseCST) = sigma_theo;
    
    % empty containers
    pC_sim_all = nan(nGaborCST, 1);
    pC_theo_all = pC_sim_all;
    mu_sim_N_all = pC_sim_all;
    mu_sim_G_all = pC_sim_all;
    
    for iGaborCST = 1:nGaborCST
        
        gaborCST_ln = gaborCST_ln_all(iGaborCST);
        
        % create Gabor (and prewhiten it)
%         gabor = exp_CreateGabor(stim, gaborCST_ln);
        gabor = template*gaborCST_ln; % to reduce computational need
        gabor_prw = fxn_prewhiten(gabor, filter_fft);
        
        % calculate the predicted center of gabor-present trials (S+N)
        mu_theo_G = nansum(gabor(:) .*template(:)); % use white noise, hence no need to prewhitened
        mu_theo_G_prw = dot(gabor_prw, template_prw); mu_theo_G_prw = real(mu_theo_G_prw) + imag(mu_theo_G_prw);
        
        if mode_prw
            mu_theo_G = mu_theo_G_prw;
        end
        criterion_opt = mu_theo_G/2; % the optimal criterion (half of the distance)
        
        % 2. calulate pC of simulated trials
        RV_stimG_allT = nan(1, nTrials);
        RV_stimN_allT = RV_stimG_allT;
        correctness_allT = RV_stimG_allT;
        iTgt_allT = RV_stimG_allT;
        resp_allT = RV_stimG_allT;
        
        for iTrial = 1:nTrials
            % create noise
            noise = randn(aper_psz, aper_psz); % white noise (not scaled by noise contrast yet)
            noise_band = ifft2(filter_fft.*fft2(noise), 'symmetric').*noiseCST; % bandpass noise
            noise = noise.*noiseCST; % white noise (scaled by noise contrast)
            noise_prw = fxn_prewhiten(noise_band, filter_fft); % prewhiten the bandpass noise; but cannot be replaced by 'noise'!! since DV is calculated in the fourier domain
            
            % make stimuli
            stim_G = noise + gabor; % Gabor-present trials (S+N)
            stim_N = noise; % Gabor-absent trials (N)
            if mode_prw
                stim_G = noise_band + gabor;
                stim_N = noise_band;
            end
            stim_G_prw = noise_prw + gabor_prw;
            stim_N_prw = noise_prw;
            
            % calculate response variable(RV)
            RV_stimG = sum(stim_G.*template, 'all');
            RV_stimN = sum(stim_N.*template, 'all');
            RV_stimG_prw = dot(stim_G_prw, template_prw); RV_stimG_prw = real(RV_stimG_prw) + imag(RV_stimG_prw);
            RV_stimN_prw = dot(stim_N_prw, template_prw); RV_stimN_prw = real(RV_stimN_prw) + imag(RV_stimN_prw);
            
            % plot each trial
            if nTrials<=3
                figure('Position', [0 0 1e3 1e3])
                subplot(3, 4, 2), imagesc(gabor), colorbar, title('Gabor'), axis square, set(gca,'XTick',[]), set(gca,'YTick',[])
                subplot(3, 4, 3), imagesc(gabor_prw), colorbar, title('Prewhitened gabor')  , axis square, set(gca,'XTick',[]), set(gca,'YTick',[])
                
                subplot(3, 4, 5), imagesc(noise), colorbar, title('white noise')  , axis square, set(gca,'XTick',[]), set(gca,'YTick',[])
                subplot(3, 4, 6), imagesc(noise_band), colorbar, title('bandpass noise')  , axis square, set(gca,'XTick',[]), set(gca,'YTick',[])
                subplot(3, 4, 7), imagesc(noise_prw), colorbar, title('Prewhitened bandpass noise')  , axis square, set(gca,'XTick',[]), set(gca,'YTick',[])
                
                subplot(3, 4, 9), imagesc(stim_G), colorbar, title('Stim (noisy gabor)')  , axis square, set(gca,'XTick',[]), set(gca,'YTick',[])
                subplot(3, 4, 10), imagesc(stim_G_prw), colorbar, title('Prewhitened Stim (noisy gabor)')  , axis square, set(gca,'XTick',[]), set(gca,'YTick',[])
                subplot(3, 4, 11), imagesc(stim_N), colorbar, title('Stim (noise)')  , axis square, set(gca,'XTick',[]), set(gca,'YTick',[])
                subplot(3, 4, 12), imagesc(stim_N_prw), colorbar, title('Prewhitened Stim (noise)')  , axis square, set(gca,'XTick',[]), set(gca,'YTick',[])
            end
            
            if mode_prw
                RV_stimG = RV_stimG_prw;
                RV_stimN = RV_stimN_prw;
            end
            
            RV_stimG_allT(iTrial) = RV_stimG;
            RV_stimN_allT(iTrial) = RV_stimN;
            
            % decide whether gabor is prs
            if rand>.5, iTgt=0; RV_stim = RV_stimN;
            else, iTgt=1; RV_stim = RV_stimG;
            end
            
            % decide response
            resp = RV_stim > criterion_opt;
            if iTgt, if resp==1, correctness=1; else, correctness=0; end
            else, if resp==1, correctness=0; else, correctness=1; end
            end
            iTgt_allT(iTrial) = iTgt;
            resp_allT(iTrial) = resp;
            correctness_allT(iTrial) = correctness;
            
        end %iTrial
        
        %% histogram of the simulated trials
        if nGaborCST<=10
            figure('Position', [0 0 600 300]), hold on
            histogram(RV_stimN_allT, 'FaceColor', 'b','FaceAlpha',.3, 'Normalization', 'probability')
            histogram(RV_stimG_allT, 'FaceColor', 'r','FaceAlpha',.3, 'Normalization', 'probability')
            xline(median(RV_stimN_allT), 'b', 'linewidth', 2);
            xline(median(RV_stimG_allT), 'r', 'linewidth', 2);
            xline(mu_theo_N, 'b--', 'linewidth', 2);
            xline(mu_theo_G, 'r--', 'linewidth', 2);
            
            xline(criterion_opt, 'k-', 'linewidth', 2);
            %             xlim([-20, 60]), ylim([0, .2])
            legend({'N', 'S+N', 'Median_{N}', 'Median_{S+N}','mu_{N}', 'mu_{S+N}', 'Opt criterion'})
            title(sprintf('Gabor cst=%.0f%%, Noise cst=%.0f%%, sigma_{theo}=%.2f\nPRS trials: %.2f +- %.2f // ABS trials: %.2f +- %.2f\npC = %.2f%% %.2f%%', ...
                gaborCST_ln*100, noiseCST*100, sigma_theo, ...
                median(RV_stimN_allT), std(RV_stimN_allT), mean(RV_stimG_allT), std(RV_stimG_allT), ...
                mean(correctness_allT)*100))
        end
        pC_sim_all(iGaborCST) = mean(correctness_allT);
        pC_theo_all(iGaborCST) = normcdf(criterion_opt, mu_theo_N, sigma_theo); % pC in gabor-pr/abs trials should be the same, hence, here only consider gabor-abs trials
        
    end % iGaborCST
    
    %% estimate contrast threshold
    [pC_min, iopt] = min(abs(pC_sim_all - perf_thresh));
    thresh_ln_allN(iNoiseCST) = gaborCST_ln_all(iopt);
    
    % plot simulated and predicted pC as a function of Gabor contrast
    if nNoiseCST<=10
        figure, hold on
        plot(gaborCST_log_all, pC_sim_all, '-bo')
        plot(gaborCST_log_all, pC_theo_all, 'k-', 'LineWidth', 2)
        
        xlabel('Gabor contrast (log)')
        ylabel('Accuracy')
        ylim([.45, 1])
        yline(.5, 'k-');
        yline(perf_thresh, 'r-');
        legend({'Simulated pC', 'Theoretical pC'}, 'Location', 'best')
        title(sprintf('Noise SD = %.2f (ntrials=%d)', noiseCST, nTrials))
        
        set(findall(gcf, '-property', 'LineWidth'), 'LineWidth',1.5)
        set(findall(gcf, '-property', 'fontsize'), 'fontsize',20)
    end
    
end % iNoise

%%
function [x_prw_fft] = fxn_prewhiten(x, filter_fft)
% prwhiten any image (x) and the filter mask in the fourier domain (filter_fft)

nHalf = size(x,1)/2;

x_fft = fft2(x); % convert the image to the fourier domain
x_fft = x_fft(1:nHalf, :); % only look at two of the four quandrants
x_fft_vec = x_fft(:); % convert to a vector

filter_fft = filter_fft(1:nHalf, :); % only look at two of the four quandrants
filter_fft_vec = filter_fft(:); % convert to a vector
ind = filter_fft_vec>0; % select the non-zero part of the mask
filter_fft_vec_s = filter_fft_vec(ind);
x_fft_vec_s = x_fft_vec(ind);

x_prw_fft = x_fft_vec_s./filter_fft_vec_s; % whiten the image

% figure
% subplot(1, 3, 1), imagesc(real(x_fft)), colorbar, axis square
% subplot(1, 3, 2), imagesc(filter_fft), colorbar, axis square
% subplot(1, 3, 3), imagesc(filter_fft.*real(x_fft)), colorbar, axis square

end

function gabor = exp_CreateGabor(stim, cst_ln)

[x1,y1] = meshgrid(1:stim.aper_psz, 1:stim.aper_psz); % size is the number of pixels
x2 = x1-mean(x1(:));
y2 = y1-mean(y1(:));
x3 = Scale(x2)*stim.aper_sz-stim.aper_sz/2; % Scale: perform an affine scaling to put data in range [0-1].
y3 = Scale(y2)*stim.aper_sz-stim.aper_sz/2; % gabor_sz is sz*ppd

ori = stim.targetOri / 180 * pi;
nx = x3*cos(ori) + y3*sin(ori);

if ~isfield(stim, 'phase'), phase =  rand * 2 *pi; else, phase = stim.phase; end

carrier = cos(nx * stim.gaborSF * 2 * pi + phase); % gaborSF_ppd is 0.0625
modulator = exp(-((x3 / stim.gaborSD).^2)-((y3 / stim.gaborSD).^2)); % gaborSD is 25.6
gabor = carrier .* modulator * cst_ln;
end

function output = Scale(input)
% Perform an affine scaling to put data in range [0-1].
minval = min(input(:));
maxval = max(input(:));
output = (input - minval) ./ (maxval-minval);
end
