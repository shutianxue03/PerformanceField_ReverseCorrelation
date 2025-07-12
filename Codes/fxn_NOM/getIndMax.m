function [maxIV, maxORI, maxSF_log] = getIndMax(IV2D_ORI_nonselected, ORI_bound, axis_tuning, plotFlag)
IV2D = IV2D_ORI_nonselected(ORI_bound(1):ORI_bound(2), :);
IV2D = IV2D_ORI_nonselected;

maxIV = max(IV2D(:));

[imaxORI, imaxSF] = find(IV2D_ORI_nonselected == maxIV);
maxORI = axis_tuning{1}(imaxORI);
maxSF_log = axis_tuning{2}(imaxSF);

% force to be one value
maxORI=maxORI(1);
maxSF_log = maxSF_log(1);

if plotFlag
    hold on
    imagesc(IV2D_ORI_nonselected), plot(imaxSF, imaxORI, 'r*')
    yline(ORI_bound(1), 'r-', 'linewidth', 2);
    yline(ORI_bound(2), 'r-', 'linewidth', 2);
    colorbar
end