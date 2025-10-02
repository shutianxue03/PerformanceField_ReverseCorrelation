

load('data.mat')


flag_GLM = 'multivariate';
[kernel2D_mul] = SX_RC(filtersSF_all, filtersOri_all, e3D, resp, flag_GLM);

flag_GLM = 'univariate';
[kernel2D_uni] = SX_RC(filtersSF_all, filtersOri_all, e3D, resp, flag_GLM);

%% Plot the kernels
axis_tuning{1} = filtersOri_all - 90;
axis_tuning{2} = log2(filtersSF_all);
axisTicks_tuning = {-90:45:90, linspace(axis_tuning{2}(1), axis_tuning{2}(end), 5)}; % ticks (SF is on log scale)
axisTL_tuning = {axisTicks_tuning{1}, round(2.^axisTicks_tuning{2}, 2)}; % label (SF is on linear scale)

figure('Position', [100 100 1200 600])
subplot(1,2,1)
imagesc(axis_tuning{2}, axis_tuning{1}, kernel2D_mul), axis square
xline(1, 'r-'); % log gabor SF
yline(0, 'r-'); % Gabor ori
xlabel('Spatial Frequency'), ylabel('Orientation')
xticks(axisTicks_tuning{2}), xticklabels(axisTL_tuning{2})
yticks(axisTicks_tuning{1}), yticklabels(axisTL_tuning{1})
title('Multivariate RC Kernel')
colorbar

subplot(1,2,2)
imagesc(axis_tuning{2}, axis_tuning{1}, kernel2D_uni), axis square
xline(1, 'r-'); % log gabor SF
yline(0, 'r-'); % Gabor ori
xlabel('Spatial Frequency'), ylabel('Orientation')
xticks(axisTicks_tuning{2}), xticklabels(axisTL_tuning{2})
yticks(axisTicks_tuning{1}), yticklabels(axisTL_tuning{1})
title('Univariate RC Kernel')
colorbar


set(findall(gcf, '-property', 'fontsize'), 'fontsize', 15)
set(findall(gcf, '-property', 'linewidth'), 'linewidth',2)  

sgtitle('RC kernel derived using different GLMs'); 


