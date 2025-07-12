clear all, clc, close all

addpath(genpath('fxn_exp'))

mode_prw = 1; % 1=using bandpass noise and prewhiten all stim
nTrials = 1e2;
nGaborCST = 5; gaborCST_log_all = linspace(-3, 0, nGaborCST); gaborCST_ln_all = 10.^gaborCST_log_all;
nNoiseCST = 2; noiseCST_all = linspace(.2, 1, nNoiseCST);
cst_ln_template = 100/100;
gaborSF = 2; SF_low = 1; SF_high = 4;
gaborSD = .8;

ppd = 32;
aper_sz = 3;
aper_psz = aper_sz*ppd;
nyquist  = ppd/2;
targetOri = 90;

stim.aper_sz=aper_sz;
stim.aper_psz=aper_psz;
stim.targetOri=targetOri;
stim.phase=0;
stim.gaborSF=gaborSF;
stim.gaborSD = gaborSD;
stim.gabor_sz = aper_sz;

perf_thresh = .7; perfThresh_all = perf_thresh;

visual.bgColor = .5;
visual.ppd = ppd;

% bandpass filter (Fourier domain)
freq     = linspace(0, nyquist, floor(aper_psz/2)+1);
faxis    = [freq, fliplr(freq(2:ceil(aper_psz/2))) ];
[fgrid_x, fgrid_y] = meshgrid(faxis, faxis);
[~, fgrid]    = cart2pol(fgrid_x, fgrid_y);
filter_fft = double(fgrid >= SF_low & fgrid <= SF_high);

% create & prewhiten the template (a full-contrast Gabor)
phase_template = 0;
template = exp_CreateGabor(stim, cst_ln_template);
template_prw = fxn_prewhiten(template, filter_fft);
limit0to1 = @(x) min(max(x,0),1);


%%
sigma_sim_allN = nan(1, nNoiseCST);
sigma_theo_allN = sigma_sim_allN;

for iNoiseCST = 1:nNoiseCST
    
    noiseCST = noiseCST_all(iNoiseCST);
    fprintf('Noise#%d/%d (%.3f)...\n', iNoiseCST, nNoiseCST, noiseCST)
    
    % theoretical variance
    if mode_prw
        %         sigma_theo = sqrt(sumsqr(template_prw))*noiseCST;
%         sigma_theo = sum(abs(template_prw))*noiseCST;
%         sigma_theo = (sum(real(template_prw)) + sum(imag(template_prw)))*noiseCST;
        sigma_theo = sum(real(template_prw)) + sum(imag(template_prw))*noiseCST;
        sigma_theo = sqrt(sum(abs(template_prw)))*noiseCST;
        
%         sigma_theo = sqrt(sum(real(template_prw).^2+imag(template_prw).^2))*noiseCST;
%         sigma_theo = sqrt(sumsqr(real(template_prw)))*noiseCST;
    else
        sigma_theo = sqrt(sum(template.^2, 'all'))*noiseCST;
    end
    sigma_theo_allN(iNoiseCST) = sigma_theo;
    
    mu_theo_N = 0;
    
    pC_sim_all = nan(nGaborCST, 1);
    mu_sim_N_all = pC_sim_all;
    mu_sim_G_all = pC_sim_all;
    
    for iGaborCST = 1:nGaborCST
        
        gaborCST_ln = gaborCST_ln_all(iGaborCST);
        
        gabor = exp_CreateGabor(stim, gaborCST_ln);
        mu_theo_G = nansum(gabor(:) .*template(:));
        gabor_prw = fxn_prewhiten(gabor, filter_fft);
        % why does fft of a Gabor has such low value (max=1) than noise (max=40)
%         mu_theo_G_prw = sum(abs(gabor_prw(:) .*template_prw(:)));
        mu_theo_G_prw = sum(real(gabor_prw.*template_prw)) + sum(imag(gabor_prw.*template_prw));
        %         mu_theo_G_prw = dot(real(gabor_prw), real(template_prw)) + dot(imag(gabor_prw), imag(template_prw)); % the same as sum(abs(complex1.*complex2))
        
        if mode_prw
            mu_theo_G = mu_theo_G_prw;
        end
        criterion_opt = mu_theo_G/2; % optimal criterion
        
        % 2. calulate pC of simulated trials
        DV_stimG_allT = nan(1, nTrials);
        DV_stimN_allT = DV_stimG_allT;
        correctness_allT = DV_stimG_allT;
        iTgt_allT = DV_stimG_allT;
        resp_allT = DV_stimG_allT;
        
        for iTrial = 1:nTrials
            % create fxn_prewhitened bandpass noise
            noise = randn(aper_psz, aper_psz);
            %             noise_band = ifftshift(ifft2(filter_fft.*fftshift(fft2(noise)), 'symmetric'));
            noise_band = ifft2(filter_fft.*fft2(noise), 'symmetric').*noiseCST;
            noise = noise.*noiseCST;
            noise_prw = fxn_prewhiten(noise_band, filter_fft); % cannot be replaced by 'noise'!! since DV is calculated in the fourier domain
            
            % make texture
            stim_G = noise + gabor;
            stim_N = noise;
            if mode_prw
                stim_G = noise_band + gabor;
                stim_N = noise_band;
            end
            stim_G_prw = noise_prw + gabor_prw;
            stim_N_prw = noise_prw;
            
            % calculate internal variable
            DV_stimG = sum(stim_G.*template, 'all');
            DV_stimN = sum(stim_N.*template, 'all');
            
%             DV_stimG_prw = sum(abs(stim_G_prw.*template_prw));
%             DV_stimN_prw = sum(abs(stim_N_prw.*template_prw));
            
%             DV_stimG_prw = (sum(real(stim_G_prw)) + sum(imag(template_prw)));
%             DV_stimN_prw = (sum(real(stim_N_prw)) + sum(imag(template_prw)));

            DV_stimG_prw = sum(real(stim_G_prw.*template_prw)) + sum(imag(stim_G_prw.*template_prw));
            DV_stimN_prw = sum(real(stim_N_prw.*template_prw)) + sum(imag(stim_N_prw.*template_prw));
            
            %             DV_stimG_prw = dot(real(stim_G_prw), real(template_prw)) + dot(imag(stim_G_prw), imag(template_prw));
            %             DV_stimN_prw = dot(real(stim_N_prw), real(template_prw)) + dot(imag(stim_N_prw), imag(template_prw));
            
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
                DV_stimG = DV_stimG_prw;
                DV_stimN = DV_stimN_prw;
            end
            
            DV_stimG_allT(iTrial) = DV_stimG;
            DV_stimN_allT(iTrial) = DV_stimN;
            
            % decide whether gabor is prs
            if rand>.5 % will be replaced by preset trial sequence
                iTgt=0; DV_stim = DV_stimN;
            else, iTgt=1; DV_stim = DV_stimG;
            end
            
            % decide response
            resp = DV_stim > criterion_opt;
            if iTgt
                if resp==1, correctness=1; else, correctness=0; end
            else
                if resp==1, correctness=0; else, correctness=1; end
            end
            iTgt_allT(iTrial) = iTgt;
            resp_allT(iTrial) = resp;
            correctness_allT(iTrial) = correctness;
            
        end %iTrial
        
        %% plot all trials of one signal cst
        if nGaborCST<=10
            figure('Position', [0 0 600 300]), hold on
            histogram(DV_stimN_allT, 'FaceColor', 'b','FaceAlpha',.3, 'Normalization', 'probability')
            histogram(DV_stimG_allT, 'FaceColor', 'r','FaceAlpha',.3, 'Normalization', 'probability')
            xline(median(DV_stimN_allT), 'b', 'linewidth', 2);
            xline(median(DV_stimG_allT), 'r', 'linewidth', 2);
            xline(median(mu_theo_N), 'b--', 'linewidth', 2);
            xline(median(mu_theo_G), 'r--', 'linewidth', 2);
            
            xline(criterion_opt, 'k-', 'linewidth', 2);
            %             xlim([-20, 60])
            ylim([0, .2])
            title(sprintf('Gabor cst=%.0f%%, Noise cst=%.0f%%, sigma_{theo}=%.2f\nPRS trials: %.2f +- %.2f // ABS trials: %.2f +- %.2f\npC = %.2f%% %.2f%%', ...
                gaborCST_ln*100, noiseCST*100, sigma_theo, ...
                median(DV_stimN_allT), std(DV_stimN_allT), mean(DV_stimG_allT), std(DV_stimG_allT), ...
                mean(correctness_allT)*100))
        end
        mu_sim_N_all(iGaborCST) = median(DV_stimN_allT);
        mu_sim_G_all(iGaborCST) = median(DV_stimG_allT);
        sigma_sim_N_all(iGaborCST) = std(DV_stimN_allT);
        sigma_sim_G_all(iGaborCST) = std(DV_stimG_allT);
        pC_sim_all(iGaborCST) = mean(correctness_allT);
        pC_theo_all(iGaborCST) = normcdf(criterion_opt, mu_theo_N, sigma_theo); % pC in gabor-pr/abs trials should be the same, hence, here only consider gabor-abs trials
        
    end % iGaborCST
    
    %% estimate thresh
    %     [pC_min, iopt] = min(abs(pC_sim_all - perf_thresh));
    [pC_min, iopt] = min(abs(pC_theo_all - perf_thresh));
    thresh_ln_allN(iNoiseCST) = gaborCST_ln_all(iopt);
    
    % plot
    if nNoiseCST<=10
        figure, hold on
        plot(gaborCST_log_all, pC_sim_all, '-bo')
        plot(gaborCST_log_all, pC_theo_all, 'k-', 'LineWidth', 2)
        
        ylim([.45, 1])
        yline(.5, 'k-');
        yline(perf_thresh, 'r-');
        legend({'Simulated pC', 'Theoretical pC'})
        title(sprintf('Noise SD = %.2f\n PRESS ANY KEY TO CONTINUE', noiseCST))
        
    end
    
end % iNoise

thresh_ln_allN*100
% noise CST = .2,   .4, .6, .8, 1
% thresh [prw]: 13.8950   37.2759   37.2759   61.0540   61%
% thresh [non-prw]: 0.9541    1.9307    2.9471    3.9069    4.4984%

%%
function [x_prw_fft, x_prw] = fxn_prewhiten(x, filter_fft)
nQuadrant = size(x,1)/2;
% convert the input patch to fft first, then divide it by the filter-mask (only a quandrant and non-zero part)
x_fft = fft2(x);
% x_fft = x_fft(1:nQuadrant, 1:nQuadrant);
x_fft = x_fft(1:nQuadrant, :);
x_fft_vec = x_fft(:);

% filter_fft = filter_fft(1:nQuadrant, 1:nQuadrant);
filter_fft = filter_fft(1:nQuadrant, :);
filter_fft_vec = filter_fft(:);
ind = filter_fft_vec>0;
filter_fft_vec_s = filter_fft_vec(ind);
x_fft_vec_s = x_fft_vec(ind);

x_prw_fft = x_fft_vec_s./filter_fft_vec_s;

% figure
% subplot(1, 3, 1), imagesc(real(x_fft)), colorbar, axis square
% subplot(1, 3, 2), imagesc(filter_fft), colorbar, axis square
% subplot(1, 3, 3), imagesc(filter_fft.*real(x_fft)), colorbar, axis square

end
