 
% check whether two ways of averaging the 2D matrix over subjects are the same
% way 1: average first, then take marg
% way 2: take marg first, then averaged over subj

n = 100;
nsubj = 10;

data2D = randn(n, n, nsubj);
% data2D2 = randn(n);

% marginalized
data2D_marg = squeeze(mean(data2D,1)); % marginalized over teh 1st dim

% 2D matrix averaged over subj
data2D_mean = squeeze(mean(data2D,3));

% average first, then take marg
data2D_mean_marg = mean(data2D_mean,1); % marginalized over teh 1st dim

% take marg first, then averaged over subj
data_marg_mean = mean(data2D_marg, 2)';

figure, hold on
plot(1:n, data2D_mean_marg, 'k-', 1:n, data_marg_mean, 'r--')
corr([data2D_mean_marg', data_marg_mean'])
% The same thing!



