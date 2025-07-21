% SX_RC7_fitting.m

% Last updated by Shutian Xue on 07/21/2025
%
% Description:
%   This function fits tuning functions to marginalized orientation or spatial frequency kernels.
%   It uses constrained optimization (fmincon with MultiStart) to estimate model parameters for each feature and type
%   The function returns fitted curves, predicted kernels, estimated parameters, and R2 metrics for each line and type.
%
% Inputs:
%   ifeature         - Feature index (1=ORI, 2=SF)
%   marg             - Marginal tuning data (nLines x nTypes x nfilters)
%   ub_full          - Upper bounds for parameters
%   lb_full          - Lower bounds for parameters
%   ifamily_perF     - Model family indices per feature
%   paramInd_perF    - Parameter index combinations per feature
%   flag_standEnergy - Flag for energy normalization
%
% Outputs:
%   marg2      - Interpolated/fitted marginal data
%   margPred2  - Model predictions for marginal data
%   margParams - Estimated model parameters
%   margR2     - R2 goodness-

function [marg2, margPred2, margParams, margR2] = SX_RC7_fitting(ifeature, marg, ub_full, lb_full, ifamily_perF, paramInd_perF, flag_standEnergy)
% load('params_RC.mat') % load some common parameters

global flag_cutMapping flag_interpolate

%-------------%
SX_RC1_setting
%-------------%

% Define fitting settings
options = optimoptions('fmincon','MaxIterations',5000,'Display','off');
nRep = 20; % number of reps to find the global minimum
sampParam0 = @(lb, ub, n) rand(1,n) .* (ub-lb) + lb; % random sample starting points between lb and ub

% Obtain parameter indices and names
nLines = size(marg, 1);% number of multiples of interpolated pints
[nLines, nTypes, nfilters] = size(marg);
if ndims(marg)==2, [nLines, nfilters] = size(marg); nTypes = 1; end

% Obtain family and parameter indices
ifamily = ifamily_perF(ifeature);
paramInd = paramInd_perF{ifeature};
nameModelParams = namesParams_all{ifamily};
nparams_full = length(nameModelParams);

% Obtain the x-axis
if ifeature==1
    axis_ln = axis_tuning{ifeature};
else
    axis_ln = 2.^axis_tuning{ifeature};
end

% Interpolate the x-axis if needed
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

% Preallocate variables
marg2 = nan(nLines, nTypes, nfilters);
margPred2 = nan(nLines, nTypes, nfilters);
margParams = nan(nLines, nTypes, nparams_full);
margR2 = nan(nLines, nTypes);

% Loop through each type of data (1=signal-present, 2=noise-present, 3=both)
for iType = 1:nTypes

    % Loop through each line
    for iLine = 1:nLines
        % Interpolate marginalized kernels
        marg_noninterp = marg;
        marg_ = squeeze(marg(iLine, iType, :));

        % Select one line
        if ndims(marg)==2, marg_ = marg(iLine, :); end
    
            if flag_interpolate>1
            marg_ = interp1(axis_ln_noninterp, squeeze(marg(iLine, iType, :))', axis_ln, 'spline');
        end
        marg2_(iLine, :) = marg_;
    end

    % Prepare bounds and initial parameters
    ub = [];
    lb = [];
    params0 = [];

    % Organize ub, lb and param0 given a certain param combination
    for iparam = 1:nparams_full
        if paramInd(iparam) == 1
            ub = [ub, ones(1, nLines)*ub_full(iparam)];
            lb = [lb, ones(1,nLines)*lb_full(iparam)];
            params0 = [params0, sampParam0(lb_full(iparam), ub_full(iparam), nLines)]; % set param0
        else
            ub = [ub, ub_full(iparam)];
            lb = [lb, lb_full(iparam)];
            params0 = [params0, sampParam0(lb_full(iparam), ub_full(iparam), 1)]; % set param0
        end
    end

    % Define the optimization problem
    fitMode = 1; % 1=ML, 2=MAP, 3=MAP with energy normalization
    fxn_minimizeDev_ML = @(kernelParams) MC_getDev_IDVD(paramInd, kernelParams, axis_ln, marg2_, ifamily, fitMode); % output of MC_getDev_IDVD is deviance
    problem_ML = createOptimProblem('fmincon','objective', fxn_minimizeDev_ML,'x0',params0,'lb',lb,'ub',ub,'options',options);
    ms_ML = MultiStart('StartPointsToRun', 'bounds','UseParallel', 1, 'Display', 'off');
    params_est = run(ms_ML, problem_ML, nRep);
    
    % Reform estimated params from a vector to 2 x nparams_full format
    params_est_ = nan(nLines, nparams_full);
    ii=1;
    for iparam = 1:nparams_full
        params_est_(1,iparam) = params_est(ii);
        params_est_(2,iparam) = params_est(ii+paramInd(iparam));
        ii = ii+nLines^paramInd(iparam);
    end

    % Predict 2D kernel sensitivity
    margPred2_ = MC_predKernel(paramInd, axis_ln, 2, nparams_full, params_est, ifamily);

    % Compute R2
    for iLine = 1:nLines
        margR2_(iLine) = 1-sumsqr(marg2_(iLine,:)-margPred2_(iLine,:))/sumsqr(marg2_(iLine,:)-mean(marg2_(iLine,:)));
    end

    %% Store results
    marg2(:, iType, :) = marg2_;
    margParams(:, iType, :) = params_est_;
    margPred2(:, iType, :) = margPred2_;
    margR2(:, iType) = margR2_;

end % end of iType





