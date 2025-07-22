function RCplot_energy(filterResp, itgt, freqStep, plotMode)

% plotMode = 1 (by SF filter), = 2 (by trial#)

if length(size(filterResp)) >2
    nAllTrials = size(filterResp,3);
    energy_SFMarg = squeeze(mean(filterResp, 2));
else
    nAllTrials = size(filterResp,1);
    energy_SFMarg = filterResp;
end

ntrials = nAllTrials/2;
energy_prs = energy_SFMarg(boolean(itgt'), :);
energy_abs = energy_SFMarg(boolean(1-itgt'), :);

switch plotMode
    case 1 % x-axis is SF channels
        plot(log(freqStep), energy_prs, 'color', [1,0,0, .1], 'handlevisibility', 'off')
        plot(log(freqStep), mean(energy_prs), 'color', [.75,0,0], 'linewidth', 2)
        plot(log(freqStep), energy_abs, 'color', [.5, .5, .5, .1], 'handlevisibility', 'off')
        plot(log(freqStep), mean(energy_abs), 'k', 'linewidth', 2)
        xlabel('SF channel')
        xticks(log([1,2,4])), xticklabels([1,2,4])
        ylim([-4,4])
        
    case 2 % x-axis is # trial
        plot(1:ntrials, energy_prs, 'color', [1,0,0,.1])
        plot(ntrials+1:nAllTrials, energy_abs, 'color', [.5, .5, .5, 1])
        plot(1:ntrials, mean(energy_prs, 2), 'k', 'linewidth', 2)
        plot(ntrials+1:nAllTrials, mean(energy_abs, 2), 'k', 'linewidth', 2)
        xticks([1,ntrials,nAllTrials]), xticklabels([1,ntrials,nAllTrials])
        xlim([1, nAllTrials])
        xlabel('trial #')
        xticks([1, nAllTrials/2, nAllTrials]), xticklabels([1, nAllTrials/2, nAllTrials])
end



