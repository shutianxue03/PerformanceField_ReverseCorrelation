function [nrmse_w_sd, nrmse_w_range] = fxn_getNRMSE(nData, pred, data)
% fxn_getNRMSE
% Weighted NRMSE per iteration.
% Inputs:
%   nData : [nCond x 1] number of trials per condition (or weights proxy)
%   pred : [nRep x nCond] predictions
%   data : [nRep x nCond] data
%
% Outputs:
%   nrmse_w_sd    : [nRep x 1] weighted RMSE normalized by weighted SD of d
%   nrmse_w_range : [nRep x 1] weighted RMSE normalized by range of d

% ---------------- checks ----------------
assert(isvector(nData), 'n must be a vector [nCond x 1].');
nData = nData(:);

assert(ismatrix(pred) && ismatrix(data), 'p and d must be 2D arrays [nRep x nCond].');
assert(all(size(pred) == size(data)), 'p and d must have the same size.');
assert(size(pred,2) == numel(nData), 'numel(n) must match size(p,2).');

nRep = size(data,1);

% ---------------- weights ----------------
% Use sqrt(n) weights and normalize to mean=1 (your convention)
w = sqrt(nData);
w = w / mean(w, 'omitnan');     % scalar normalization
W = repmat(w(:)', nRep, 1); % [nRep x nCond]

% ---------------- weighted RMSE ----------------
err2 = (data - pred).^2;

% Zero-out weights where err2 is NaN (so denom matches effective data)
W_eff = W;
W_eff(isnan(err2)) = 0;

num_rmse = sum(W_eff .* err2, 2);     % [nRep x 1]
den_rmse = sum(W_eff, 2);            % [nRep x 1]

rmse_w = sqrt(num_rmse ./ den_rmse);  % correct weighted RMSE
rmse_w(den_rmse == 0) = NaN;

% ---------------- (2A) Normalize by weighted SD of data per iteration ----------------
% Also handle NaNs consistently: weights zero where d is NaN
W_d = W;
W_d(isnan(data)) = 0;

ybar_w = sum(W_d .* data, 2) ./ sum(W_d, 2);     % [nRep x 1]
ybar_w(sum(W_d,2) == 0) = NaN;

sd_w = sqrt( sum(W_d .* (data - ybar_w).^2, 2) ./ sum(W_d, 2) ); % [nRep x 1]
sd_w(sum(W_d,2) == 0) = NaN;

% guard against sd==0
nrmse_w_sd = rmse_w ./ max(sd_w, eps);

% ---------------- (2B) Normalize by range of data per iteration ----------------
rangeData = max(data, [], 2, 'omitnan') - min(data, [], 2, 'omitnan'); % [nRep x 1]

% guard against range==0
nrmse_w_range = rmse_w ./ max(rangeData, eps);

end