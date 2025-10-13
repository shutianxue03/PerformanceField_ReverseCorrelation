
% plot if the IV distribution (raw and the transformed version) passed the
% normality test (i.e., swtest)

% h_allSubj: nsubj x nLoc8 x nModelA x nModelB x ni x 4 x 2
%% plot for each subj
iModelA=1;
iModelB=1;
for ii = 1:ni
    for isubj = 1:nsubj
        
        for iLoc = 1:nLoc8
            hh = squeeze(h_allSubj(isubj, iLoc, iModelA, iModelB, ii, :, :));
            
        end
        
    end
end
%%

figure('Position', [0 0 1600 600])

iplot = [4:7, 11:14, 18:21];
for isubj = 1:nsubj
    hh = squeeze(mean(h_allSubj(isubj, :, :, :, :, :, :)==0, 1));
    hh_allSubj(isubj, :, :) = hh;
    subplot(3,7,iplot(isubj))
    bar(hh)
    xticks(1:length(namesNormalityTest))
    xticklabels(namesNormalityTest), xtickangle(45)
    ylim([0,1])
    title(subjList{isubj})
end

%%
subplot(3,7,[1:3, 8:10, 15:17]), hold on
bar(squeeze(mean(hh_allSubj,1)))
buffer=.15;
errorbar((1:nNormTests)-buffer, squeeze(mean(hh_allSubj(:,:, 1))), squeeze(std(hh_allSubj(:,:, 1)))/sqrt(nsubj), 'k.', 'CapSize', 0)
errorbar((1:nNormTests)+buffer, squeeze(mean(hh_allSubj(:,:, 2))), squeeze(std(hh_allSubj(:,:, 2)))/sqrt(nsubj), 'k.', 'CapSize', 0)

xticks(1:length(namesNormalityTest))
xticklabels(namesNormalityTest), xtickangle(45)
ylim([0,1])
ylabel('Prop. of passing swtest')
legend({'PRS', 'ABS'})
title(sprintf('n = %d ave.', nsubj))

set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
