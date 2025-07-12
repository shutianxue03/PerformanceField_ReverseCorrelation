
% check the energy profile of a gabor of diff SFs

nfiltersOri = length(filtersOri_all);
nfiltersSF = length(filtersSF_all);


%%
for gaborSF = 1:4
    % gaborSF=2;
    stim.gaborSF=gaborSF;
    gabor = exp_CreateGabor(stim, 1, 0);
    noiseP = exp_CreateFilteredNoise(noise);
    pp = gabor;
    
    energy = nan(nfiltersOri, nfiltersSF);
    stimPhase = energy;
    
    for ifilterOri = 1:nfiltersOri
        for ifilterSF = 1:nfiltersSF
            [a,b] = SX_sim04_computeEnergy(pp , pp(1,1), filter_sin{ifilterOri, ifilterSF}, filter_cos{ifilterOri, ifilterSF}, 0, 1, filtersSF_all(ifilterSF), filtersOri_all(ifilterOri));
            energy(ifilterOri, ifilterSF) = a;
            stimPhase(ifilterOri, ifilterSF) = b;
        end
    end
    
    %% stimulus phase and 2D energy
    figure('Position', [0 0 1000 400])
    subplot(1,3,1), imshow(pp)
    
    subplot(1,3,2), hold on
    imagesc(filtersOri_all-90, filtersSF_all_log, stimPhase')
    colorbar
    % caxis([-3, 3])
    xline(0, 'r-', 'linewidth', 2);
    yline(1, 'r-', 'linewidth', 2);
    xticks(-80:40:80)
    xlim([-80,80])
    yticks([0,1,2]), yticklabels([1,2,4])
    ylabel('SF (cpd)'), xlabel('ORI (deg)')
    axis square
    title('Stim Phase')
    
    subplot(1,3,3), hold on
    imagesc(filtersOri_all-90, filtersSF_all_log, energy')
    colorbar
    xline(0, 'r-', 'linewidth', 2);
    yline(1, 'r-', 'linewidth', 2);
    xticks(-80:40:80)
    xlim([-80,80])
    yticks([0,1,2]), yticklabels([1,2,4])
    ylabel('SF (cpd)'), xlabel('ORI (deg)')
    axis square
    title('Energy (unnormalized)')
    
    sgtitle(['Gabor SF = ', num2str(gaborSF)])
    
    set(findall(gcf, '-property', 'FontSize'), 'FontSize',18)
    
    %%
    energy_norm = (energy-mean(energy(:)))/std(energy(:));
    
    figure
    
    subplot(2,2,1), hold on
    plot(filtersSF_all_log, mean(energy))
    xline(1, 'r');
    xticks([0,1,2]), xticklabels([1,2,4])
    xlabel('SF (cpd)'), ylabel('Energy')
    ylim([0,1])
    
    subplot(2,2,2), hold on
    plot(filtersOri_all-90, mean(energy,2))
    xline(0, 'r');
    xticks(-80:40:80)
    xlabel('ORi (deg)'), ylabel('Energy')
    ylim([0,1])
    
    subplot(2,2,3), hold on
    plot(filtersSF_all_log, mean(energy_norm))
    xline(1, 'r');
    yline(0, 'r');
    ylim([-1, 1])
    xticks([0,1,2]), xticklabels([1,2,4])
    xlabel('SF (cpd)'), ylabel('Normalized energy')
    
    subplot(2,2,4), hold on
    plot(filtersOri_all-90, mean(energy_norm,2))
    xline(0, 'r');
    yline(0, 'r');
    ylim([-1, 1])
    xticks(-80:40:80)
    xlabel('ORi (deg)'), ylabel('Normalized energy')
    
    sgtitle(['Gabor SF = ', num2str(gaborSF)])
    set(findall(gcf, '-property', 'FontSize'), 'FontSize',18)
    
end