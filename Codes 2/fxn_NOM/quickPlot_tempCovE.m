

load('Data_OOD/params_RC')

[nORI, nSF] = size(tempCovE);
figure('Position', [3e3 0 1e3 400])
subplot(1,3,1), hold on
RCplot_2Dkernel(template', [])
maxIV = max(template(:));
[imaxORI_, imaxSF_] = find(template == maxIV);
maxORI_ = axis_tuning{1}(imaxORI_);
maxSF_log_ = axis_tuning{2}(imaxSF_);
xline(maxORI_, 'm--', 'linewidth', 2);
yline(maxSF_log_, 'm--', 'linewidth', 2);
title('TEMPATE')

subplot(1,3,2), hold on
RCplot_2Dkernel(e', [])
maxIV = max(e(:));
[imaxORI_, imaxSF_] = find(e == maxIV);
maxORI_ = axis_tuning{1}(imaxORI_);
maxSF_log_ = axis_tuning{2}(imaxSF_);
xline(maxORI_, 'm--', 'linewidth', 2);
yline(maxSF_log_, 'm--', 'linewidth', 2);
title('Energy profile')

subplot(1,3,3), hold on
RCplot_2Dkernel(tempCovE', [])
xline(maxORI, 'm--', 'linewidth', 2);
yline(maxSF_log, 'm--', 'linewidth', 2);
switch convolveType
    case 1, title('Dot product')
    case 2 , title('Convolution')
end


set(findall(gcf, '-property', 'FontSize'), 'FontSize',18)
sgtitle(sprintf('Trial#%d, ORI max = %d, SF max = %.2f', itrial, maxORI, round(2.^maxSF_log, 2)))


