
% function MCplot_allMode
for MCmode = 1:3 % 1 or 2 or 3
    
    % general params
    SX_RC1_setting
    
    itype = 3; % ONLY look at kernels derived from ALL trials; 1=PRS; 2=ABS; 3=ALL
    iLoc_all = [1,8]; nLoc2 = length(iLoc_all);
    subjList = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA', 'CS', 'DT'};
    nblocks_allSubj = [200, 240, 210, 240, 220, 220, 215, 220, 205, 210, 205, 205];
    nsubj = length(subjList);
    
    %% model setting
    if MCmode == 3, nn = nIC; else, nn = 1; end
    % names of params
    namesParams = namesParams_all{ifamily};
    nparams_full = length(namesParams);
    
    % decide number of models & the param index (1 = specified, 0 = shared)
    if sum(ifamily == [3, 7]) % the truncation term is always left free
        nmodels = 2^(nparams_full-1);
        paramInd_all_ = fxn_getParamInd(nparams_full-1);
        paramInd_all = [paramInd_all_, ones(nmodels,1)];
    else
        nmodels = 2^nparams_full;
        paramInd_all = fxn_getParamInd(nparams_full);
    end
    
    % extract ub and lb
    ub_full = ub_full_all{ifamily};
    lb_full = lb_full_all{ifamily};
    
    %% extract x (y is extracted for each subj)
    if (ifeature == 2) && (~sum(ifamily == [4, 5, 9])) % fitting SF except (4) raised gaussian (5) double exponential (9) gaussian to SF
        axis_tuning_ln = 2.^axis_tuning{ifeature};
    else, axis_tuning_ln = axis_tuning{ifeature};
    end
    nfilters = length(axis_tuning_ln);
    
    %% load MC files
    nameFileMC = sprintf('Data_MC/n%d_%s_Family%d_mode%d.mat', nsubj, namesFeature{ifeature}, ifamily, MCmode);
    load(nameFileMC)
    fprintf('Loaded\n')
    
    %% plot the dev/IC & freq averaged across subj
    MCplot_perMCmode
    
    %% plot the data and prediction
    MCplot_pred
    
end
