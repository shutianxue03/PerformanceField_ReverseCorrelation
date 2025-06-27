%% look at Yuna's data
% gpeak is the peak CS 
% fpeak is the peak SF. 
% They are both in scale of log10. 
% Each row is each participant, and the columns are in the order of upper / lower / horizontal.

clear all, clc
load data_neutral_forSX.mat
% nsubj x nLoc
[nsubj, nLoc] = size(fpeak_neut);

indLoc = repmat((1:nLoc)', [nsubj, 1]);

data = fpeak_neut; % peak SF
data = gpeak_neut; % peak CS

text_ANOVA = print_nANOVA({'Loc'}, data(:), {indLoc(:)}, nsubj);

HM = data(:,3);
VM = mean(data(:, [1,2]), 2);

[~, p, ~, stats] =ttest(HM, VM);
text_HVA = sprintf('t=%.2f, p=%.4f (%d/%d)', stats.tstat, p, sum(HM>VM), nsubj);
[~, p, ~, stats] =ttest(data(:, 2), data(:, 1));
text_VMA= sprintf('t=%.2f, p=%.4f (%d/%d)', stats.tstat, p, sum(data(:, 2)>data(:, 1)), nsubj);

fprintf('%s\nHVA: %s\nVMA: %s\n', text_ANOVA, text_HVA, text_VMA)

