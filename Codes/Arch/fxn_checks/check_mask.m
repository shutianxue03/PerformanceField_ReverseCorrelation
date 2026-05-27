
% check different parameters when creating mask

sinsz_ppd_all = (.2:.2:1)*ppd;
sinpower_all = [0,3,5,7];
iplot=1;

figure('Position', [0 0 1e3, 1e3])
for sinsz_ppd = sinsz_ppd_all
    for sinpower = sinpower_all 
        stim.sin_psz = sinsz_ppd;
        stim.sinpower = sinpower;
        
        [patch, full_objsiz] = exp_CreateCircularApertureSin(stim);
        subplot(length(sinsz_ppd_all), length(sinpower_all),iplot), imagesc(patch), axis square
        colorbar
        title(sprintf('sin sz = %.1f\nsin power = %d', round(sinsz_ppd/ppd,1), sinpower))
        iplot = iplot+1;
    end
end