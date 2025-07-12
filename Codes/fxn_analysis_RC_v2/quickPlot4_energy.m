e2D = energy2D_allT_perComb; title_ = 'Unnormalized';
% e2D = energy2D_norm_allT_perComb; title_ = 'Normalized';

%% marginalized energy
for iLoc = 1%:nLoc8
    plotEnergy1D(e2D{iLoc})
    sgtitle(sprintf('%s energy -  %s', title_, namesLocComb{iLoc}))
end

%% stim phase
figure('Position', [0 200 1500 600])
for iLoc = 1:nLoc8
    subplot(2,4,iLoc)
    imagesc(squeeze(mean(stimPhase2D_allT_perComb{iLoc},1)))
    axis square
    colorbar
    title(namesLocComb{iLoc})
end
sgtitle('Stim phase')
set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)

%%
function plotEnergy1D(e2D)
% for a full version, see quickPlot_energy_arch

figure('Position', [0 500 1000 800])
nAllTrials = size(e2D, 1);
ntrials = nAllTrials/2;
% e2D_ave = squeeze(mean(e2D, 1));
% e2D_std = squeeze(std(e2D, 1));

for ifeature = 1:2
    if ifeature == 1, e_marg = squeeze(mean(e2D, 3));
    else, e_marg = squeeze(mean(e2D, 2));
    end
    
    % energy as a fxn of trials
    subplot(2,2,1+(ifeature-1)*2), hold on
    plot(1:ntrials, e_marg(1:ntrials, :), 'color',[1,0,0,.1]),
    plot(ntrials+1:ntrials*2, e_marg(ntrials+1:end, :), 'color',[0,0,1,.1])
    xlabel('left is PRS right is ABS')
    
    % energy as a fxn of ORI/SF
    subplot(2,2,2+(ifeature-1)*2), hold on
    plot(e_marg(1:ntrials, :).', 'color',[1,0,0,.1]),
    plot(e_marg(ntrials+1:end, :).' , 'color',[0,0,1,.1]),
end
set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)



end