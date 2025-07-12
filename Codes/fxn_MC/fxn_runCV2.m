function [param_est_allK, dev_test_allK, y_pred_allK] = fxn_runCV2(CVmode, k, ncv, yData, ifamily, options, params0, lb, ub, nrep, axis_tuning_ln)
% k:
%    LOOCV: nfilters;
%    k-fold CV: k
nfilters = length(yData);
nparams = length(params0);

% empty containers
param_est_allK = nan(k * ncv, nparams);
dev_test_allK = nan(k * ncv,1);
dev_allData_allK = dev_test_allK;
y_pred_allK = nan(k * ncv, nfilters);

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
        
        % x
        x_test = axis_tuning_ln(ind_test);
        x_train = axis_tuning_ln(ind_train);
        % y
        y_test = yData(ind_test);
        y_train= yData(ind_train);
        
        % MCmode serves as fitMode (1) for RSS (2) for nLL
        fitMode = 2;
        %%
        [params_est, dev_test, dev_allData, y_test_pred, y_data_pred] = fxn_MC_getEST(yData, x_test, y_test, x_train, y_train, ifamily, fitMode, options, params0, lb, ub, nrep, axis_tuning_ln);
        
%         if rand<.1, figure, hold on, plot(x_train, y_train, 'ko'), plot(axis_tuning_ln, y_data_pred, 'k-'), plot(x_test, y_test_pred, '*r'), end % quick plot
        param_est_allK(ik_cv, :) = params_est;
        dev_test_allK(ik_cv) = dev_test;
        dev_allData_allK(ik_cv) = dev_allData;
        y_pred_allK(ik_cv, :) = y_data_pred;
    end % ik
end % end of icv

%% helper fxn
    function [params_est, dev_test, dev_allData, y_test_pred, y_data_pred] = fxn_MC_getEST(yData, x_test, y_test, x_train, y_train, ifamily, fitMode, options, params0, lb, ub, nrep, axis_tuning_ln)
        fxn_getDev = @(kernelParams) MC_getDev_CV2(kernelParams, x_train, y_train, ifamily, fitMode);
        problem = createOptimProblem('fmincon','objective', fxn_getDev,'x0',params0,'lb',lb,'ub',ub,'options',options);
        ms = MultiStart('StartPointsToRun', 'bounds', 'UseParallel', 1, 'Display', 'off');
        params_est = run(ms, problem, nrep);
        % validate
        dev_test = MC_getDev_CV2(params_est, x_test, y_test, ifamily, fitMode);
        dev_allData = MC_getDev_CV2(params_est, axis_tuning_ln, yData, ifamily, fitMode);
        % predict
        y_test_pred = predSFkernel(x_test, ifamily, params_est, 0);
        y_data_pred = predSFkernel(axis_tuning_ln, ifamily, params_est, 0);
    end
end
