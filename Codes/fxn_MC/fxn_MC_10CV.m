function fxn_MC_10CV
% repeating each CV
param_est_allK = nan(k * ncv, nparamsAll);
dev_test_allK = nan(k * ncv,1);
dev_allData_allK = dev_test_allK;

for icv = 1:ncv
    for ik = 1:k % in each iteration, params are estimated based on the training group and validation based on the test group
        ik_cv = k*(icv-1)+ik;
        if ik==1, istart = 1; else, istart = sz_cumsum(ik-1)+1;end
        iend = sz_cumsum(ik);
        ntest = length(istart:iend);
        ntrain = nfilters - ntest;
        % empty containers
        x_test = nan(2, ntest);
        y_test  = x_test;
        x_train = nan(2, ntrain);
        y_train = x_train;
        for iline = 1:2
            % get the index for the test/train group
            ind_test = ind(istart:iend, iline);
            ind_train = 1:nfilters; ind_train(ind_test) = [];
            assert(length(ind_test)+length(ind_train) == nfilters)
            % x
            x_test(iline, :) = axis_tuning_ln(ind_test);
            x_train(iline, :) = axis_tuning_ln(ind_train);
            % y
            y_test(iline, :) = yData(iline, ind_test);
            y_train(iline, :) = yData(iline, ind_train);
        end % end of iline
        
        % MCmode serves as fitMode (1) for RSS (2) for nLL
        fitMode = MCmode;
        [params_est, dev_test, dev_allData, y_test_pred, y_data_pred] = fxn_MC_getEST(paramInd, nparams_full, yData, x_test, y_test, x_train, y_train, ifamily, fitMode, options, params0, lb, ub, nrep, axis_tuning_ln);
        param_est_allK(ik_cv, :) = params_est;
        dev_test_allK(ik_cv) = dev_test;
        dev_allData_allK(ik_cv) = dev_allData;
    end % ik
end % end of icv
param_est_all{imodel} = mean(param_est_allK, 1);
dev_ave(imodel) = mean(dev_test_allK);
dev_std(imodel)  = std(dev_test_allK);
end