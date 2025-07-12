

function [marg2, margPred2, margParams, margR2] = SX_RC7_fitting(ifeature, marg, ub_full, lb_full, ifamily_perF, paramInd_perF, flag_standEnergy)
% load('params_RC.mat') % load some common parameters

global flag_cutMapping flag_interpolate
SX_RC1_setting
options = optimoptions('fmincon','MaxIterations',5000,'Display','off');
nrep = 20; % number of reps to find the global minimum
sampParam0 = @(lb, ub, n) rand(1,n) .* (ub-lb) + lb;
nlines = size(marg, 1);% number of multiples of interpolated pints
[nlines, ntypes, nfilters] = size(marg);
if ndims(marg)==2, [nlines, nfilters] = size(marg); ntypes = 1; end

ifamily = ifamily_perF(ifeature);
paramInd = paramInd_perF{ifeature};
nameModelParams = namesParams_all{ifamily};
nparams_full = length(nameModelParams);

% get x
if ifeature==1
    axis_ln = axis_tuning{ifeature};
else
    axis_ln = 2.^axis_tuning{ifeature};
end
% axis_ln_noninterp = axis_ln;
% if flag_interpolate>1
%     if ifeature==1 % only interpolate the center
%         nNotItp = sum(axis_ln<=-30);
%         axis_ln_itp = linspace(axis_ln(nNotItp), axis_ln(end-nNotItp+1), flag_interpolate); 
%         axis_ln = [axis_ln(1:nNotItp-1), axis_ln_itp, axis_ln(end-nNotItp+2:end)];
%     else
%         axis_ln = linspace(axis_ln(1), axis_ln(end), flag_interpolate); 
%     end
% end
nfilters = length(axis_ln);


% empty containers
marg2 = nan(nlines, ntypes, nfilters);
margPred2 = nan(nlines, ntypes, nfilters);
margParams = nan(nlines, ntypes, nparams_full);
margR2 = nan(nlines, ntypes);

for itype = 1:ntypes

    for iline = 1:nlines
        % interpolate marg
        marg_noninterp = marg;
        marg_ = squeeze(marg(iline, itype, :));
        if ndims(marg)==2, marg_ = marg(iline, :); end
        if flag_interpolate>1
            marg_ = interp1(axis_ln_noninterp, squeeze(marg(iline, itype, :))', axis_ln, 'spline');
        end
        marg2_(iline, :) = marg_;
    end
    
    ub = [];
    lb = [];
    params0 = [];
    % organize ub, lb and param0 given a certain param combination
    rng default % For reproducibility
    for iparam = 1:nparams_full
        if paramInd(iparam) == 1
            ub = [ub, ones(1, nlines)*ub_full(iparam)];
            lb = [lb, ones(1,nlines)*lb_full(iparam)];
            params0 = [params0, sampParam0(lb_full(iparam), ub_full(iparam), nlines)]; % set param0
        else
            ub = [ub, ub_full(iparam)];
            lb = [lb, lb_full(iparam)];
            params0 = [params0, sampParam0(lb_full(iparam), ub_full(iparam), 1)]; % set param0
        end
    end
    
    fitMode = 1;
    fxn_minimizeDev_ML = @(kernelParams) MC_getDev_IDVD(paramInd, kernelParams, axis_ln, marg2_, ifamily, fitMode); % output of MC_getDev_IDVD is deviance
    problem_ML = createOptimProblem('fmincon','objective', fxn_minimizeDev_ML,'x0',params0,'lb',lb,'ub',ub,'options',options);
    ms_ML = MultiStart('StartPointsToRun', 'bounds','UseParallel', 1, 'Display', 'off');
    params_est = run(ms_ML, problem_ML, nrep);
    % reform params to 2 x nparams_full
    params_est_ = nan(nlines, nparams_full);
    ii=1;
    for iparam = 1:nparams_full
        params_est_(1,iparam) = params_est(ii);
        params_est_(2,iparam) = params_est(ii+paramInd(iparam));
        ii = ii+nlines^paramInd(iparam);
    end
    
    % predict
    margPred2_ = MC_predKernel(paramInd, axis_ln, 2, nparams_full, params_est, ifamily);
    
    % Rsquared
    for iline = 1:nlines
        margR2_(iline) = 1-sumsqr(marg2_(iline,:)-margPred2_(iline,:))/sumsqr(marg2_(iline,:)-mean(marg2_(iline,:)));
    end
    
    %% compile
    marg2(:, itype, :) = marg2_;
    margParams(:, itype, :) = params_est_;
    margPred2(:, itype, :) = margPred2_;
    margR2(:, itype) = margR2_;

end % end of itype





