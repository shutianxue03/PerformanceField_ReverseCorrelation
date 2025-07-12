
% ePRS/ABS: ntrials x nORI x nSF
% could be per loc per itgt (PRS/ABS/both)

e3D_PRS = e3D(dataMatrix(:, 5)==iLoc & dataMatrix(:, 6)==1, :, :);
e3D_ABS = e3D(dataMatrix(:, 5)==iLoc & dataMatrix(:, 6)==0, :, :);

%% non-standardized
figure('Position', [0 300 2e3 300])
% mean across trials
subplot(1,5,1), imagesc(squeeze(mean(e3D_PRS,  1))), colorbar, axis square, title('[PRS] Mean across trials')
% std across trials
subplot(1,5,2), imagesc(squeeze(std(e3D_PRS,  [], 1))), colorbar, axis square, title('[PRS] SD across trials')
% mean across trials
subplot(1,5,3), imagesc(squeeze(mean(e3D_ABS,  1))), colorbar, axis square, title('[ABS] Mean across trials')
% std across trials
subplot(1,5,4), imagesc(squeeze(std(e3D_ABS,  [], 1))), colorbar, axis square, title('[ABS] SD across trials')
sgtitle(sprintf('%s-%s [non-standardized]', subjName, namesLocComb{iLoc}))

% all trials of a specific channel
iORI_signal = ceil(nORI/2); 
iSF_signal = ceil(nSF/2);

subplot(1,5,5), hold on
plot(squeeze(e3D_PRS(:, iORI_signal, iSF_signal)), 'r')
plot(squeeze(e3D_ABS(:, iORI_signal, iSF_signal)), 'b')
legend({'PRS', 'ABS'})
xlabel('Trial #')
ylabel('Energy')
title('Energy on 0 deg and 2 cpd')

%% marginalized ORI
figure('Position', [0 0 800 300])
subplot(1,2,1), hold on

iDim = 3;
plot(axis_tuning{1}, squeeze(mean(mean(e3D_PRS, 1), iDim)), 'r')
plot(axis_tuning{1}, squeeze(mean(mean(e3D_ABS, 1), iDim)), 'b')
legend({'PRS', 'ABS'})
xlabel('ORI')
xticks(axisTicks_tuning{1})
xticklabels(axisTL_tuning{1})
ylim([0, .15])
title('Marginalized ORI energy')

% marginalized SF
subplot(1,2,2), hold on

iDim = 2;
plot(axis_tuning{2}, squeeze(mean(mean(e3D_PRS, 2), 1)), 'r')
plot(axis_tuning{2}, squeeze(mean(mean(e3D_ABS, 2), 1)), 'b')
legend({'PRS', 'ABS'})
xlabel('SF')
xticks(axisTicks_tuning{2})
xticklabels(axisTL_tuning{2})
ylim([0, .15])
title('Marginalized SF energy')

sgtitle(sprintf('%s [%s]', subjName, namesLocComb{iLoc}))
set(findall(gcf, '-property', 'fontsize'), 'fontsize', 15)
set(findall(gcf, '-property', 'LineWidth'), 'LineWidth',1.5)
