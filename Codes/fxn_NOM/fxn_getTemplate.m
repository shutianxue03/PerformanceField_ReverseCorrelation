
function [template] = fxn_getTemplate(kernels2D, templateType, flag_plot)

%% get templates
template_raw = kernels2D;
% [nSF, nORI] = size(kernels2D);
[nORI, nSF] = size(kernels2D);
t_min = min(template_raw(:));
%  reconstruct the template
template_pos = template_raw - t_min + eps; % shift to make all values positive
template_recon = mean(template_pos,1) .* mean(template_pos,2) + t_min;
% mirror the template
indMir = (nORI-1)/2;
kk_left = template_raw(1:indMir, :);
kk_right = template_raw(indMir+2:end, :);
kk_mid = template_raw(indMir+1,:);
kk_ave = (kk_left + flip(kk_right))/2;
template_mirrored = [kk_ave; kk_mid; flip(kk_ave)];
% mirror and reconstruct the template
template_rm = template_mirrored - min(template_mirrored(:)) + eps;
template_rm = mean(template_rm,1) .* mean(template_rm,2) + min(template_mirrored(:));

%%
if flag_plot
    figure('Position', [0 0 1e3 400])
    subplot(1,4,1), hold on, imagesc(template_raw), axis square, yline((nORI+1)/2); xline((nSF+1)/2); colorbar, title('raw')
    subplot(1,4,2), hold on, imagesc(template_recon), axis square, yline((nORI+1)/2); xline((nSF+1)/2);colorbar, title('reconstructed')
    subplot(1,4,3), hold on, imagesc(template_mirrored), axis square, yline((nORI+1)/2); xline((nSF+1)/2);colorbar, title('mirrored')
    subplot(1,4,4), hold on, imagesc(template_rm), axis square, yline((nORI+1)/2); xline((nSF+1)/2);colorbar, title('m&r')
end

%% set template type
switch templateType
    case 1, template = template_raw;% (1) the 2D kernel (raw)
    case 2, template = template_recon; % (2) the reconstructed 2D kernel
    case 3, template = template_mirrored; % (3) the mirrored 2D kernel
    case 4, template = template_rm; % (4) the % mirror and reconstruct the template
    case 5, template = energy2D_gabor;
end