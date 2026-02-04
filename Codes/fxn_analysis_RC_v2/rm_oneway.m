function [Fvalue, eta2_partial] = rm_oneway(data)
% Repeated-measures one-way ANOVA F for Condition (within-subject)
% Standard partition: SS_cond, SS_error (subject x condition interaction)

[nSubj, nCond] = size(data);

grand = mean(data(:), 'omitnan');
subj_mean = mean(data, 2, 'omitnan');     % nSubj x 1
cond_mean = mean(data, 1, 'omitnan');     % 1 x nCond

SS_cond = nSubj * sum((cond_mean - grand).^2, 'omitnan');
SS_subj = nCond * sum((subj_mean - grand).^2, 'omitnan');

% Total SS
SS_tot  = sum((data - grand).^2, 'all', 'omitnan');

% Error term in RM one-way: subject x condition interaction
SS_err = SS_tot - SS_cond - SS_subj;

df_cond = nCond - 1;
df_err  = (nSubj - 1) * (nCond - 1);

MS_cond = SS_cond / df_cond;
MS_err  = SS_err  / df_err;

Fvalue = MS_cond / MS_err;

% Compute the effect size: partial eta-squared
eta2_partial = SS_cond/(SS_cond+SS_err);
end