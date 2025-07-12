
%% extract ub and lb
ub_full = ub_full_all{ifamily};
lb_full = lb_full_all{ifamily};

%% extract x (y is extracted for each subj)
if (ifeature == 2) && (~sum(ifamily == [4, 5, 9])) % fitting SF except (4) raised gaussian (5) double exponential (9) gaussian to SF
    axis_tuning_ln = 2.^axis_tuning{ifeature};
else, axis_tuning_ln = axis_tuning{ifeature};
end
nfilters = length(axis_tuning_ln);

%% Fig 1. plot the dev/IC & freq averaged across subj
if flagPlotComp
    %-------------------%
    MCplot1_perMCmode
    %-------------------%
end

%% Fig 2. plot the data and prediction
% flag_plot=0;
%-------------------%
MCplot2_pred
%-------------------%

