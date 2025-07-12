% function quickPlot_energy

namesFeature = {'ORI', 'SF'};
namesType = {'PRS', 'ABS', 'BOTH'};
namesLoc2D = {'Fovea', 'Left', 'Upper', 'Right', 'Lower'};
namesLocComb = {'Fovea', 'LHM', 'UVM', 'RHM', 'LVM', 'HM', 'VM', 'Peri'};
namesSDT = {'d''', 'criterion', 'RT (sec)'};
namesEnergy = {'WHOLE', 'NOISE', 'NOISE-masked'};
namesNorm = {'Unnorm', 'Norm'};
namesRespType = {'Hit rate', 'FA rate', 'pC'};
namesTitles1 = {'Fovea', 'Horizontal', 'Lower', 'Left'}; % before 'vs.', color-coded
namesTitles2 = {'Periphery', 'Vertical', 'Upper', 'Right'}; % after 'vs.', color-coded
namesTitles3 = {'', 'Meridian', 'Vertical Meridian', 'HM'}; % in black
namesCI = {'mean', 'mean norm', 'var', 'var norm'};
%%
ntrials = size(e2D, 1)/2;

%% stimulus phase
% if normFLAG == 0
%     figure
%     imagesc(squeeze(mean(stimPhase2D_allT,1)).')
%     colorbar
%     xline(0, 'r-', 'linewidth', 2);
%     yline(1, 'r-', 'linewidth', 2);
% %     xticks(-80:40:80)
% %     xlim([-80,80])
% %     yticks([0,1,2]), yticklabels([1,2,4])
% %     ylabel(namesFeature{2}), xlabel(namesFeature{1})
%     axis square
% %     title('Stim Phase')
% %     title(['Stim Phase-',namesLoc2D{iLoc}, '-',namesEnergy{energyFlag}], 'FontSize',25)
%     set(findall(gcf, '-property', 'FontSize'), 'FontSize',18)
% end


%% 2D plot
% namesMeanSD = {'Mean', 'SD'};
% figure('Position', [0 0 600 1000])
% 
% for ii = 1:2 % mean and std
%     if ii==1 % mean (ideally=0)
%         caxisMin = -.1; caxisMax = .1;
%     else % std (ideally=1)
%         caxisMin = .9; caxisMax = 1.1;
%     end
%     
%     for itype = 1:ntypes % 1=PRS, 2=ABS, 3=BOTH
%         switch itype
%             case 1, itrial_start = 1; itrial_end = ntrials;
%             case 2, itrial_start = ntrials+1; itrial_end = 2*ntrials;
%             case 3, itrial_start = 1; itrial_end = 2*ntrials;
%         end
%         
%         if ii==1, e2D_ = squeeze(mean(e2D(itrial_start:itrial_end, :, :),1)); % nOri x nSF
%         else, e2D_ = squeeze(std(e2D(itrial_start:itrial_end, :, :),[], 1)); % nOri x nSF
%         end
%         
%         subplot(ntypes, 2, ii+(itype-1)*2)
%         imagesc(filtersOri_all-90, filtersSF_all_log, e2D_')
%         colorbar
%         xline(0, 'r-', 'linewidth', 2);
%         yline(1, 'r-', 'linewidth', 2);
%         xticks(-80:40:80)
%         xlim([-80,80])
%         yticks([0,1,2]), yticklabels([1,2,4])
%     
%         ylabel(namesFeature{2}), xlabel(namesFeature{1})
%         axis square
%         
%         if normFLAG, caxis([caxisMin, caxisMax]), end
%         
%         title(['Energy-',namesType{itype}, '-', namesMeanSD{ii}])
%     end
% end
% set(findall(gcf, '-property', 'FontSize'), 'FontSize',18)
% sgtitle([namesLoc2D{iLoc}, '-', namesEnergy{energyFlag}], 'FontSize',25)

%% 1D
figure('Position', [0 500 1000 800])
e2D_ave = squeeze(mean(e2D, 1));
e2D_std = squeeze(std(e2D, 1));

for ifeature = 1:2
    if ifeature == 1
%         filters_all = filtersOri_all-90;
        e_marg = squeeze(mean(e2D, 3));
        %         e2D_ave_marg = mean(e2D_ave,2);
        %         e2D_std_marg = mean(e2D_std,2);
    else
%         filters_all = filtersSF_all;
        e_marg = squeeze(mean(e2D, 2));
        %         e2D_ave_marg = mean(e2D_ave,1);
        %         e2D_std_marg = mean(e2D_std,1);
    end
%     nfilters = length(filters_all);
%     xticks_ = [1,ceil(nfilters/2), nfilters];
    
    % energy as a fxn of trials
    subplot(2,2,1+(ifeature-1)*2), hold on
    plot(1:ntrials, e_marg(1:ntrials, :), 'color',[1,0,0,.1]),
    plot(ntrials+1:ntrials*2, e_marg(ntrials+1:end, :), 'color',[0,0,1,.1])
    xlabel('left is PRS right is ABS')
    
    % energy as a fxn of ORI/SF
    subplot(2,2,2+(ifeature-1)*2), hold on
    plot(e_marg(1:ntrials, :).', 'color',[1,0,0,.1]),
    plot(e_marg(ntrials+1:end, :).' , 'color',[0,0,1,.1]),
%     xticks(xticks_), xticklabels(filters_all(xticks_))
%     xlabel(namesFeature{ifeature})
    
    % mean and SD
    %     subplot(2,2,3+(ifilter-1)*3), hold on
    %     plot(filters_all, e2D_ave_marg, 'ro')
    %     plot(filters_all, e2D_std_marg, 'r+')
    %     plot(mean(e_marg(ntrials+1:end, :)), 'bo', 'MarkerSize', 10)
    %     plot(std(e_marg(ntrials+1:end, :)), 'b+', 'MarkerSize', 10)
    %     xticks(xticks_), xticklabels(filters_all(xticks_))
    %     legend('prs mean', 'prs sd', 'abs mean', 'abs sd', 'Location', 'southwest')
    
    set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
end

% sgtitle([namesNorm{normFLAG+1}, '-', namesLocComb{iLoc}, '-',namesEnergy{energyFlag}], 'FontSize',25)
