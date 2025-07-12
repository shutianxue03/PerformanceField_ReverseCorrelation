

figure('Position', [3e3 500 1000 500])
yfit_odd = norminv(yfit_);

%% data
subplot(2,2,1), hold on
% data
plot(x(resp==1), resp(resp==1), 'co')
plot(x(resp==0), resp(resp==0), 'mo')
% data correctlresp categorized
indCorr_resp1 = resp(yfit_>=.5)==1;
indCorr_resp0 = resp(yfit_<=.5)==0;
plot(x(indCorr_resp1), resp(indCorr_resp1), 'k+', 'MarkerSize', 1)
plot(x(indCorr_resp0), resp(indCorr_resp0), 'k+', 'MarkerSize', 1)

% xline(0);
% xlim([-4, 4])
% ylim([-.2,1.2]), yticks(0:.5:1)
xlabel('Standardized energy')
ylabel('resp (1=Yes 0=NO)')

%% regression line
subplot(2,2,3), hold on
plot(x, yfit_, 'k.')
yline(0);
% xlim([-4, 4])
ylabel('p(resp=YES)')

%% histogram
subplot(2,2,[2,4]), hold on
x_resp1 = x(boolean(resp));
x_resp2 = x(boolean(1-resp));
histogram(x_resp1, 'FaceColor', 'c', 'EdgeColor', 'none', 'FaceAlpha', .3)
histogram(x_resp2, 'FaceColor', 'm', 'EdgeColor', 'none', 'FaceAlpha', .3)
xline(mean(x_resp1), 'c-');
xline(mean(x_resp2), 'm-');
xline(median(x_resp1), 'c--');
xline(median(x_resp2), 'm--');
% xlim([-3, 3])
legend({'Resp=1', 'Resp=0'})
xlabel('Standardized energy')
ylabel('Count')

%%
sgtitle(sprintf('Ori = %d, SF = %.2f\nslope = %.2f, intercept = %.2f\nTjur R^2 = %.4f, p = %.2f, pCat = %.2f', ...
    filtersOri_all(iORI)-90, filtersSF_all(iSF), ...
    beta(2), beta(1), ...
    R2_Tjur_, stats.p(2), pCat_))
set(findall(gcf, '-property', 'FontSize'), 'FontSize',18)
