
% plotFlag = 0;

%% way 1: eq 6 in AE2006
[template, ~] = exp_CreateGabor(stim, stim.gaborCST, 0);

DV_all = nan(1, ntrialsAll); % decision var

for ss = 1:ntrialsAll
    DV_all(ss) = sum(patch_both{ss}(:) .* template(:)) + randn * noise.noiseProp;
end

DV_all = (DV_all - mean(DV_all))/std(DV_all);

figure, hold on
histogram(DV_all(1:ntrials),'FaceColor', 'r', 'normalization', 'probability')
histogram(DV_all(ntrials+1:end), 'FaceColor', 'b', 'normalization', 'probability')
ylabel('decision var')
legend('PRS', 'ABS')
xline(criterion_true, 'k', 'linewidth', 2);

behav = DV_all > criterion_true;

%% way 2: RC
% lumiBG = template(1,1);
% energy_temp = nan(1,nfilters);
% for ff = 1:nfilters
%     [a,b] = SX_sim04_computeEnergy(template , lumiBG, SFfilter_sin{ff}, SFfilter_cos{ff}, 0, ss<=ntrials, filterSF_all(ff));
%     energy_temp(ff) = a;
% end
% energy_temp_norm = energy_temp - mean(energy_temp);
% gaussNoise = randn(1,ntrials*2) * noise.noiseProp;
% 
% resp_allTrials = energy_temp_norm * energy_norm';
% resp_noisy = resp_allTrials + gaussNoise;

%% way 3: original one

% SX_sim05_getBehavior
% gaussNoise = randn(1,ntrials*2) * noise.noiseProp;
% 
% resp_allTrials = energy_temp_norm * energy_norm';
% resp_noisy = resp_allTrials + gaussNoise;


%%
% if plotFlag
%     figure('Position',[0 200 400 1200])
%     suptitle('Responses')
%     
%     subplot(4,1,1), hold on
%     plot(filtersSF_all, resp')
%     plot(filtersSF_all([1,end]), [0 0], 'Color', [.5 .5 .5])
%     xlabel('SF channel (cpd)'), title('kernel x energy')
%     
%     subplot(4,1,2), hold on
%     plot(filterSF_all, mean(resp(1:ntrials, :)))
%     plot(filterSF_all, mean(resp(ntrials+1:end, :)))
%     plot(filterSF_all([1,end]), [0 0], 'Color', [.5 .5 .5])
%     xlabel('SF channel (cpd)'), title('kernel x energy (averaged)')
%     legend('tgt-prs', 'tgt-abs')
%     
%     subplot(4,1,3), hold on
%     stem(resp_noisy')
%     plot([0,ntrials*2], [0 0], 'Color', [.5 .5 .5])
%     title('noise added, not normalized')
%     
%     subplot(4,1,4), hold on
%     stem(resp_norm')
%     ylimit = ylim; ymax = ylimit(2);ymin = ylimit(1);
%     plot([ntrials, ntrials], [ymin, ymax], 'k-')
%     plot([0,ntrials*2], [0 0], 'Color', [.5 .5 .5])
%     title(sprintf('noise added, normalized (criterion = %d)', stim.criterion))
%     xlabel('trial #')
%     text(ntrials, 0.15, 'tgt-prs trials        tgt-absent trials','HorizontalAlignment', 'center')
%     
%     set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
%     
% end
