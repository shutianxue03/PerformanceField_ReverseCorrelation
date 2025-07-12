function [ave, lb, ub, SEM_neg, SEM_pos] = getCI(mat, errType, dim, flagSqueeze, ntrialsProp, CI_level)

% [ave, lb, ub, SEM_neg, SEM_pos] = getCI(mat, errType, dim, flagSqueeze, ntrialsProp, CI_level)
%   errType: (1) get 68% confidence interval (default) (2) get SEM
%   dim: the dimension to be averaged (default = 1)
%   flagSqueeze: whether squeeze the output (default = 1, squeeze)
%   ntrialsProp: the prop of ntrials per subj (default = 1, do NOT weight average based on trial prop)
%   CI_level: default = 0.68

%% define default values
if nargin < 2, errType = 1; end
if nargin < 3, dim = 1; end
if nargin < 4, flagSqueeze = 1; end
if nargin < 5, ntrialsProp = 1; end
if nargin < 6, CI_level = .68; end

%%
if errType == 2, nsubj = size(mat, 1); end

switch errType
    case 1 % get median and 68% CI
        ave = nanmedian(mat, dim);
        lb = quantile(mat, .5-CI_level/2, dim);
        ub = quantile(mat, .5+CI_level/2, dim);
        SEM_neg = ave - lb;
        SEM_pos = ub - ave;
        
    case 2 % get mean and SEM
        ave = squeeze(nanmean(mat, dim));
        SEM_neg = squeeze(nanstd(mat, [], dim))/sqrt(nsubj);
        SEM_pos = SEM_neg;
        lb = ave - SEM_neg;
        ub = ave + SEM_pos;
end

if flagSqueeze
    ave = squeeze(ave);
    SEM_neg = squeeze(SEM_neg);
    SEM_pos = squeeze(SEM_pos);
    lb = squeeze(squeeze(lb));
    ub = squeeze(ub);
end

