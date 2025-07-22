

SF_lb = 1.22;
SF_lb_log = log2(SF_lb);
SF_ub = 3.28;
SF_ub_log = log2(SF_ub);
SFend = 20; %round(stim.aper_sz*4); % the 
xaxis = (1:SFend)/stim.aper_sz;
xaxis_log = log2(xaxis);

figure('Position', [0 200 1500 400])
for tt = 1:2 % prs and abs
    
    for iLoc = 1:nLoc8
        fft1D_mean = squeeze(fft1D_mean_both(iLoc, tt, :));
        if nsubj==1, fft1D_std = zeros(length(fft1D_mean),1);end
        
        subplot(2, nLoc8, nLoc8*(tt-1)+iLoc), hold on, box on
        stem(xaxis_log, fft1D_mean(2:SFend+1))
        errorbar(xaxis_log, fft1D_mean(2:SFend+1), fft1D_std(2:SFend+1), 'k', 'linestyle', 'none')
        xline(log2(2), 'r');
        xlim([-1.1, 3])
        poly = polyshape([SF_lb_log, SF_ub_log, SF_ub_log, SF_lb_log], [0,0,max(fft1D_mean), max(fft1D_mean)]);
        plot(poly, 'FaceColor','red','FaceAlpha',0.1, 'Edgecolor', 'w')
        xticks(log2(1:4)), xticklabels(1:4)
        xlabel('SF (cpd)')
%         ylim([0,20])
        if tt==1, title(namesLocComb{iLoc}), end
    end
end

set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
sgtitle('top row PRS + bottom row ABS', 'fontsize', 25)

