function [ave, lb, ub, SEM_neg, SEM_pos] = getCI(mat, errType, dim, CI_level, flag_squeezeDims)

% [ave, lb, ub, SEM_neg, SEM_pos] = getCI(mat, errType, dim, flagSqueeze, ntrialsProp, CI_level)
%   errType: (1) get 68% confidence interval (default) (2) get SEM
%   dim: the dimension to be averaged (default = 1)
%   CI_level: default = 0.68
%   flagSqueeze: whether squeeze the output (default = 1, squeeze)

%% define default values
if nargin < 2, errType = 1; end
if nargin < 3, dim = 1; end
if nargin < 4, CI_level = .68; end
if nargin < 5, flag_squeezeDims = 1; end

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
        SEM = squeeze(nanstd(mat, [], dim))/sqrt(nsubj);
%         SEM = withinSubjErr(mat); 
        SEM_neg = SEM;
        SEM_pos = SEM;
        lb = ave - SEM_neg;
        ub = ave + SEM_pos;
end

if flag_squeezeDims
    ave = squeeze(ave);
    SEM_neg = squeeze(SEM_neg);
    SEM_pos = squeeze(SEM_pos);
    lb = squeeze(squeeze(lb));
    ub = squeeze(ub);
end

