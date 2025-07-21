% fxn_runCV.m
%
% Last updated by Shutian Xue on 07/21/2025
%
% Description:
%   This function performs cross-validation (CV) for tuning function fitting to conduct model comparison (MC) analysis.
%   It supports both k-fold and leave-one-out cross-validation (LOOCV) modes, partitions the data, 
%   fits model parameters using constrained optimization (MultiStart with fmincon), and evaluates model performance on test sets.
%   The function returns estimated parameters, test set deviance, and predicted values for each CV fold.
%
% Inputs:
%   CVmode         - Cross-validation mode (1 = k-fold, 2 = LOOCV)
%   k              - Number of folds (for k-fold CV) or number of data points (for LOOCV)
%   ncv            - Number of cross-validation repetitions
%   nfilters       - Number of data points (filters)
%   yData          - Observed data (2 x nfilters)
%   paramInd       - Indices of parameters to fit
%   ifamily        - Model family index
%   options        - Optimization options for fmincon
%   params0        - Initial parameter values
%   lb, ub         - Lower and upper bounds for parameters
%   nrep           - Number of MultiStart runs
%   axis_tuning_ln - Tuning axis for model prediction
%
% Outputs:
%   param_est_allK - Estimated parameters for each fold/repetition
%   dev_test_allK  - Deviance on test set for each fold/repetition
%   y_pred_allK    - Predicted values for each


function [param_est_allK, dev_test_allK, y_pred_allK] = fxn_runCV(CVmode, k, ncv, nfilters, yData, paramInd, ifamily, options, params0, lb, ub, nrep, axis_tuning_ln)
% k:
%    LOOCV: nfilters
%    k-fold CV: k

nparamsAll = sum(2.^paramInd);
nparams_full = length(paramInd);

% empty containers
param_est_allK = nan(k * ncv, nparamsAll);
dev_test_allK = nan(k * ncv,1);
dev_allData_allK = dev_test_allK;
y_pred_allK = nan(k * ncv, 2, nfilters);

for icv = 1:ncv % repetation
    if CVmode == 1 % 10-fold CV
        hpartition = cvpartition(nfilters, 'kfold', k);
        sz_perSet = hpartition.TestSize; % the number of data points in each set
        sz_cumsum = cumsum(sz_perSet); % the number of data points accumulatively
        rng('shuffle')
        ind = randperm(nfilters, nfilters);
    end
    
    for ik = 1:k % in each iteration, params are estimated based on the training group and validation based on the test group
        ik_cv = k*(icv-1)+ik;
        switch CVmode
            case 1 % k-fold CV
                % istart of the test set
                if ik==1, istart = 1; else, istart = sz_cumsum(ik-1)+1;end
                % iend of the test set
                iend = sz_cumsum(ik);
                ind_test = ind(istart:iend);
            case 2 % LOOCV
                ind_test = ik;
        end
        
        ind_train = 1:nfilters; ind_train(ind_test) = [];
        ntest = length(ind_test);
        ntrain = nfilters - ntest;
        assert(length(ind_test)+length(ind_train) == nfilters)
        
        % empty containers
        x_test = nan(2, ntest); y_test  = x_test;
        x_train = nan(2, ntrain); y_train = x_train;
        
        % could below could be put in a fxn
        for iline = 1:2
            % x
            x_test(iline, :) = axis_tuning_ln(ind_test);
            x_train(iline, :) = axis_tuning_ln(ind_train);
            % y
            y_test(iline, :) = yData(iline, ind_test);
            y_train(iline, :) = yData(iline, ind_train);
        end % end of iline
        
        % MCmode serves as fitMode (1) for RSS (2) for nLL
        fitMode = 1;
        %%
        %-------------------------------------------------%
        [params_est, dev_test, dev_allData, y_test_pred, y_data_pred] = fxn_MC_getEST(paramInd, ...
            nparams_full, yData, x_test, y_test, x_train, y_train, ifamily, fitMode, options, params0, lb, ub, nrep, axis_tuning_ln);
        %-------------------------------------------------%
        param_est_allK(ik_cv, :) = params_est;
        dev_test_allK(ik_cv) = dev_test;
        dev_allData_allK(ik_cv) = dev_allData;
        y_pred_allK(ik_cv, :, :) = y_data_pred;
    end % ik
end % end of icv

%% helper fxn
    function [params_est, dev_test, dev_allData, y_test_pred, y_data_pred] = fxn_MC_getEST(paramInd, nparams_full, yData, x_test, y_test, x_train, y_train, ifamily, fitMode, options, params0, lb, ub, nrep, axis_tuning_ln);
        %-------------------------------------------------%
        fxn_getDev = @(kernelParams) MC_getDev_CV(paramInd, kernelParams, x_train, y_train, ifamily, fitMode);
        %-------------------------------------------------%
        problem = createOptimProblem('fmincon','objective', fxn_getDev,'x0',params0,'lb',lb,'ub',ub,'options',options);
        ms = MultiStart('StartPointsToRun', 'bounds', 'UseParallel', 1, 'Display', 'off');
        [params_est, ~] = run(ms, problem, nrep);
        % validate
        dev_test = MC_getDev_CV(paramInd, params_est, x_test, y_test, ifamily, fitMode);
        dev_allData = MC_getDev_CV(paramInd, params_est, axis_tuning_ln, yData, ifamily, fitMode);
        % predict
        y_test_pred = MC_predKernel(paramInd, x_test, 2, nparams_full, params_est, ifamily);
        y_data_pred = MC_predKernel(paramInd, axis_tuning_ln, 2, nparams_full, params_est, ifamily);
    end
end
