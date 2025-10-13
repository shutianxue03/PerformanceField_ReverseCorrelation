function [R2, SSR] = getR2(y, y_pred, weights)
if isnan(weights), weights=ones(size(y)); end
weighted_mean_y = sum(weights .* y) / sum(weights);
% Weighted SSR (Sum of Squared Residuals)
SSR = sum(weights .* (y - y_pred).^2);
% Weighted SST (Total Sum of Squares)
SST = sum(weights .* (y - weighted_mean_y).^2);
% Weighted R^2
R2 = 1 - (SSR / SST);
end