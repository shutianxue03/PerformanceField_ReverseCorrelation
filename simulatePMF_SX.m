clc
% close all

SF = 2;
ORI = 1; % signal is right-tilted
perf_level = .7;

nTrials = 1e2;
cst_template = 1;
nCST = 10; cst_all = 10.^linspace(-3, 0, nCST);
nNoiseSD = 1; noiseSD_all = linspace(.5, nNoiseSD);
SF_low = 1;
SF_high = 4;

ppd = 32; % pixel per degree

gaborsz = 2.*ppd;
gaborsz = round(gaborsz/2)*2;
gaborenvelopedev = 1.*ppd;
gaborangle = [-45 45];
gaborfrequency = SF/ppd;

limit0to1 = @(x) min(max(x,0),1);

%% create the template
gaborPhase = 0; % constant phase
template_L = CreateGabor_SX(gaborsz,gaborenvelopedev,135,gaborfrequency,gaborPhase, cst_template);
template_R = CreateGabor_SX(gaborsz,gaborenvelopedev,45,gaborfrequency,gaborPhase, cst_template);
template_diff =  template_L-template_R;

%%
figure, hold on

sigma_sim_all = nan(1, nNoiseSD);
% sigma_theo_all = sigma_sim_all;
% sigma_theo_all = [nan .747 1.5171 2.2618 3.0953] .* 1e5;

f = waitbar(0,sprintf('Simulating (0/%i)...', nNoiseSD));
for iNoise = 2%2:nNoiseSD
    waitbar(iNoise/nNoiseSD, f, sprintf('Simulating (%i/%i)...', iNoise, nNoiseSD))
    % calculated the theoretical sigma
    noiseSD = noiseSD_all(iNoise);
    sigma_theo = sqrt(sumsqr(template_diff))*noiseSD; % only depends on noiseSD
    %     sigma_theo_all(iNoise) = sigma_theo;
    
    pC_sim_all = nan(1, nCST);
    pC_theo_all = nan(1, nCST);
    
    for iCST = 1:nCST
        
        cst = cst_all(iCST);
        %         cst = .051;
        %         cst = .5;
        
        % create the Gabor at the given contrast
        gabor = CreateGabor_SX(gaborsz,gaborenvelopedev,gaborangle(ORI),gaborfrequency,...
            gaborPhase, cst);
        
        %%%%%%%%%%%% simulation %%%%%%%%%%%%
        internalVar_allT = nan(1, nTrials);
        correctness_allT = internalVar_allT;
        
        for iTrial = 1:nTrials
            % create noisy Gabor
            noiseImg = randn(gaborsz, gaborsz).*noiseSD; %Make the noise and adjust its contrast
            nyquist  = ppd/2;
            freq     = linspace(0,nyquist,floor(gaborsz/2)+1);
            faxis    = [freq, fliplr(freq(2:ceil(gaborsz/2))) ];
            faxis    = fftshift(faxis);
            [fgrid_x,fgrid_y] = meshgrid(faxis,faxis);
            [~, fgrid]    = cart2pol(fgrid_x,fgrid_y);
            
            %% Draw filter for frequency domain
            filter_fft = double(fgrid > SF_low & fgrid < SF_high);
            %             filter_fft(:,1:length(freq)-1) = 0;
            
            %% Generate Gaussian noise and apply the filters
            fn    = fftshift(fft2(noiseImg));
            noiseImg = ifftshift(ifft2(fftshift(1.*filter_fft.*fn))); %changed this
            noiseImg = noiseImg/std(noiseImg(:)).*noiseSD - mean(noiseImg(:));
            
            %             noiseImg = randn(gaborsz, gaborsz) .* noiseSD;
            noise_only = noiseImg;
            noiseImg = noiseImg + gabor;
            
            % calculate the internal variable
            %             internalVar = sum(noiseImg(:).*template_diff, 'all');
            filter_fft(:,1:length(freq)-1) = 0;%fftshift(fft2(fftshift(noise_only)));
            stim_fft = fftshift(fft2(fftshift(noiseImg)));
            template_fft = fftshift(fft2(fftshift(template_diff)));
            
            flattend_nonzero = stim_fft; 
            flattend_nonzero(filter_fft~=0) = flattend_nonzero(filter_fft~=0) ./ filter_fft(filter_fft~=0);
            flattemplate_nonzero = template_fft; 
            flattemplate_nonzero(filter_fft~=0) = flattemplate_nonzero(filter_fft~=0) ./ filter_fft(filter_fft~=0);
            dv = dot(real(flattemplate_nonzero), real(flattend_nonzero)) + dot(imag(flattemplate_nonzero), imag(flattend_nonzero));
            %             internalVar = nansum((noiseImg .* template_diff) .* fftshift(ifft2(1./(((abs(fft2(fftshift(noise_only)))))))), 'all');
            internalVar = nansum(dv,'all');
            %             internalVar = nansum(noiseImg(:) .* real(fftshift(ifft2(fftshift(dv(:))))));
            
            if ORI==1 % Left
                resp = internalVar>0;
            else
                resp = internalVar<0;
            end
            if resp==1, correctness=1; else, correctness=0; end
            
            internalVar_allT(iTrial) = internalVar;
            correctness_allT(iTrial) = correctness;
        end %iTrial
        
        pC_sim_all(iCST) = mean(correctness_allT);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% calculate theoretical pC
        %         mu_theo = sum(gabor(:) .*template_diff(:));
        fft_gabor = fftshift(fft2(fftshift(gabor)));
        flattend_nonzero = fft_gabor; flattend_nonzero(filter_fft~=0) = flattend_nonzero(filter_fft~=0) ./ filter_fft(filter_fft~=0);
        dv = dot(real(flattemplate_nonzero), real(flattend_nonzero)) + dot(imag(flattemplate_nonzero), imag(flattend_nonzero));
        mu_theo = sum(dv(:));
        mu_flattemplate = mean(flattemplate_nonzero(:));
        flattemplate_quadrants = flattemplate_nonzero;
        
        %         sigma_theo = (sqrt(sumsqr(abs(flattemplate_nonzero.*.44)))).*noiseSD.* sqrt(numel(flattemplate_nonzero));
        sigma_theo = sqrt((sumsqr(real(flattemplate_nonzero(:,1:length(freq)-1))) + sumsqr(imag(flattemplate_nonzero(:,1:length(freq)-1)))).*noiseSD).*sqrt(numel(flattemplate_nonzero));
        %         sigma_theo = sigma_theo_all(iNoise);
        mu_theo_all(iCST) = mu_theo;
        sigma_theo_all(iNoise) = sigma_theo;
        if ORI==1
            pC_theo_all(iCST) = 1-normcdf(0, mu_theo, sigma_theo);
        else
            pC_theo_all(iCST) = normcdf(0, mu_theo, sigma_theo);
        end
        
    end % iCST
    
    dev =@(cst) (est_pC(10^cst, noiseSD, ppd, sigma_theo) - perf_level).^2;
    options = optimoptions('fmincon','Display','iter', 'Algorithm', 'sqp');
    thresh_est(iNoise) = (fmincon(dev, -3, [], [], [], [], -3, 0, [], options));
    
    plot(log10(cst_all), pC_sim_all, '-ro')
    plot(log10(cst_all), pC_theo_all, 'k-', 'LineWidth', 2)
    yline(perf_level, 'k-');
    if ~isnan(thresh_est(iNoise))
        xline(thresh_est(iNoise), 'k-');
    end
    
    xlabel('Log cst')
    ylabel('Accuracy')
    
    xlim([min(log10(cst_all)), max(log10(cst_all))])
    ylim([.45, 1])
    yline(.5, 'k-');
    
    pause(.5)
    
end % iNoise
close(f)

%%
function pC_theo = est_pC(cst, noiseSD, ppd, sigma_theo)

SF = 2;
ORI = 2; % signal is right-tilted
cst_template = 1;
SF_high = 4;
SF_low = 1;

%     ppd = 64; % pixel per degree
gaborsz = 2.*ppd;
gaborsz = round(gaborsz/2)*2;
gaborenvelopedev = 1.*ppd;
gaborangle = [-45 45];
gaborfrequency = SF/ppd;
gaborPhase = 0; % constant phase
template_L = CreateGabor_SX(gaborsz,gaborenvelopedev,135,gaborfrequency,gaborPhase, cst_template);
template_R = CreateGabor_SX(gaborsz,gaborenvelopedev,45,gaborfrequency,gaborPhase, cst_template);
template_diff =  template_L-template_R;

% Bandpassed noise stuff
nyquist  = ppd/2;
freq     = linspace(0,nyquist,floor(gaborsz/2)+1);
faxis    = [freq, fliplr(freq(2:ceil(gaborsz/2))) ];
faxis    = fftshift(faxis);
[fgrid_x,fgrid_y] = meshgrid(faxis,faxis);
[~, fgrid]    = cart2pol(fgrid_x,fgrid_y);
filter_fft = double(fgrid > SF_low & fgrid < SF_high);
fft_noise = filter_fft; fft_noise(:,1:length(freq)-1) = 0;
fft_template = fftshift(fft2(fftshift(template_diff)));

flattemplate_nonzero = fft_template; flattemplate_nonzero(fft_noise~=0) = flattemplate_nonzero(fft_noise~=0) ./ fft_noise(fft_noise~=0);

gabor = CreateGabor_SX(gaborsz, gaborenvelopedev, gaborangle(ORI), gaborfrequency,...
    gaborPhase, cst);
fft_gabor = fftshift(fft2(fftshift(gabor))); fft_gabor(:,1:length(freq)-1) = 0;
flattend_nonzero = fft_gabor; flattend_nonzero(filter_fft~=0) = flattend_nonzero(filter_fft~=0) ./ filter_fft(filter_fft~=0);
dv = dot(real(flattemplate_nonzero), real(flattend_nonzero)) + dot(imag(flattemplate_nonzero), imag(flattend_nonzero));
mu_theo = sum(dv(:));
%     mu_theo = sum(gabor(:) .*template_diff(:));
%     sigma_theo = sqrt(sumsqr(template_diff))*noiseSD; % only depends on noiseSD
%     sigma_theo = sqrt(sumsqr(dv(:)))*noiseSD;
if ORI==1
    pC_theo = 1-normcdf(0, mu_theo, sigma_theo);
else
    pC_theo = normcdf(0, mu_theo, sigma_theo);
    
end
end