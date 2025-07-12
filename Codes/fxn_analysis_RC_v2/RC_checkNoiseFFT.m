
% function: checkNoiseFFT
% psz = stim.aper_psz;
psz = stim.aper_psz;

dimInd = 2; % aalong which dim the 2D FFT is averaged
fft1D_mean_both = nan(nLoc8, 2, psz);

for iLoc = 1:nLoc8
    if iLoc>5, switch iLoc, case 6, iLoc_all = [2,4]; case 7, iLoc_all = [3,5]; case 8, iLoc_all = 2:5; end
    else, iLoc_all = iLoc;
    end
    ntrials_ = ntrials*length(iLoc_all);
    
    for tt = 1:2 % prs and abs
        fft2_2D = nan(ntrials_, psz, psz);
        fft2_1D = nan(ntrials_, psz);
        nn_perLoc = noise_both_perComb_p{iLoc};
        
        parfor ip = 1:ntrials_
            F = fft2(nn_perLoc{(tt-1)*ntrials_+ip});
            fft2_2D(ip, :, :) = abs(F);
            fft2_1D(ip, :) = abs(mean(F, dimInd));
        end
        
        fft2D_mean = squeeze(mean((fft2_2D), 1)); % average over all trials
        fft1D_mean = mean(fft2_1D, 1); % average over all trials
        fft1D_mean_both(iLoc, tt, :) = fft1D_mean;
    end
end



