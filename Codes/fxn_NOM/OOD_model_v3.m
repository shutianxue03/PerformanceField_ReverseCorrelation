datetime('now')
clc
% clear all % do NOT clear all, the inputs fed to the current function in the shell file will be removed
warning off
format compact

addpath(genpath('Data_OOD'))
addpath(genpath('fxn_model'))
addpath(genpath('fxn_RCplot'))
addpath(genpath('fxn_analysis_RC_v2'))
load('Data_OOD/params_RC.mat')

SX_RC1_setting

%% params to vary for each batch
isubj = 1;
iLocComb = 1;
iModelA = 1; % (1) core model (2) template is random (3) energy is random
iModelB = 2; % (1) model both additive and multiplicative internal noise (2) only additive (3) only multiplicative
minResample = 2;
flagIncludePA = 1; % whether incorporate pA in calculating nLL
flagConstrainThresh = 1;

%% general params
ni = 10;
nrep = 15; % 20
templateType = 1; % (1) raw kernel (2) positive kernel (3) reconstructed kernel
convolveType = 1; % (1) dot multiply (2) convolution
IVType = 1; % (1) sum up (2) max (3) normalization model
itype_template = 3;
nCountdown = 10;

if minResample>1, templateType=2;end

% if flagIncludePA
%     switch iModelB
%         case 1, params0 = [1, 1, 0, nan]; params_lb = [1e-5, 1e-5, 0, nan]; params_ub = [10, 10, .5, nan];
%         case 2, params0 = [1, 0, nan]; params_lb = [1e-5,  0, nan]; params_ub = [10, .5, nan];
%         case 3, params0 = [1, 0, nan]; params_lb = [1e-5,  0, nan]; params_ub = [10, .5, nan];
%     end
% else
    switch iModelB
        case 1, params0 = [1, 1, nan]; params_lb = [1e-5, 1e-5, nan]; params_ub = [10, 10, nan];
        case 2, params0 = [1, nan]; params_lb = [1e-5, nan]; params_ub = [10, nan];
        case 3, params0 = [1, nan]; params_lb = [1e-5, nan]; params_ub = [10, nan];
    end
% end


%%
SX_RC1_setting

ratio_train = 2/3; % proportion of trials in the training set (to derive 2D kernels); others are in the test set (to estimate params)
ORI_bound = [7, 23]; % -40 and 40

getAIC_ML = @(nLL, nparamsAll, ndata) 2 * nLL + 2*nparamsAll;
getAIC_LS = @(err, nparamsAll, ndata) ndata * log(err/ndata) + 2*nparamsAll;
getAICc_ML = @(nLL, nparamsAll, ndata) 2 * nLL + 2*nparamsAll + 2*nparamsAll*(nparamsAll+1)/(ndata-nparamsAll-1);
getAICc_LS = @(err, nparamsAll, ndata) ndata * log(err/ndata) + 2*nparamsAll + 2*nparamsAll*(nparamsAll+1)/(ndata-nparamsAll-1);
getBIC = @(nLL, nparamsAll, ndata) 2*nLL + nparamsAll * log(ndata);
options = optimoptions('fmincon','MaxIterations', 1e4, 'Display','off');

markers_allSubj = {'o', 's', 'd', '^','v',  '<', '>','p', 'h', '+', 'x', '-'}; % for each subj

indCountdown1 = linspace(ni/nCountdown, ni, nCountdown);
indCountdown2 = nCountdown:-1:1;

%% names
subjList_ = {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA', 'CS', 'DT'};
nblocks_allSubj_ = [200, 240, 210, 240, 220, 220, 215, 220, 205, 210, 205, 205];

namesMetrics = {'dprime',  'Criterion', 'pC', 'pHit', 'pFA', 'pA', 'pA [PRS]', 'pA [ABS]'};
nmetrics = length(namesMetrics);

namesIC = {'AIC', 'AICc', 'BIC'};
nIC = length(namesIC);

namesModelA = {'Core', 'Random template', 'Random energy'};
nModelA = length(namesModelA);

% if flagIncludePA
%     namesModelB = {'Add. + Multi. + pA buffer + thresh', 'Add. + pA buffer + thresh', 'Multi. + pA buffer + thresh'};
%     namesParamsModel_all = {{'Multi noise', 'Additive noise', 'pA buffer', 'Threshold'}, ...
%         {'Additive noise', 'pA buffer', 'Threshold'}, {'Multi noise','pA buffer', 'Threshold'}};
%     namesParams2 = {'Thresh', 'mu [PRS]', 'sigma [PRS]','mu [ABS]', 'sigma [ABS]', 'pA buffer'};
% else
    namesModelB = {'Add. + Multi. + thresh', 'Add. + thresh.', 'Multi. + thresh'};
    namesParamsModel_all = {{'Multi noise', 'Additive noise', 'Threshold'}, ...
        {'Additive noise', 'Threshold'}, {'Multi noise', 'Threshold'}};
    namesParams2 = {'Thresh', 'mu [PRS]', 'sigma [PRS]','mu [ABS]', 'sigma [ABS]'};
    
% end

nModelB = length(namesModelB);
namesParamsModel = namesParamsModel_all{iModelB};
nparams_model = length(namesParamsModel); % number of params to fit; currently: alpha and thresh
assert(nparams_model == length(params0))

nparams2 = length(namesParams2);

namesNormalityTest = {'Raw', 'Log', 'Sqrt', 'boxcox'}; % reciprocal is deleted
nNormTests = length(namesNormalityTest);
% Box, G. E. P. and Cox, D. R. (1964). An analysis of transformations

subjName = subjList_{isubj};
nblocks = nblocks_allSubj_(isubj);

%% file name & print
nameFileModelIDVD = sprintf('Data_OOD/Model_%s_Loc%d_A%dB%d_min%d_flagpA%d_thresh%d', ...
    subjName, iLocComb, iModelA, iModelB, minResample, flagIncludePA, flagConstrainThresh);

fprintf('\n%s nblocks = %d\nLoc: %s\nMODEL: %s & %s\nMin. resampling = %d\nWhether to include incorporate pA in calculating nLL: %d\nWhether constrain threshold based on criterion: %d\n', ...
    subjName, nblocks, namesLocComb{iLocComb}, namesModelA{iModelA}, namesModelB{iModelB}, minResample, ...
    flagIncludePA, flagConstrainThresh)

%% load behav data
fprintf('\nLoading the behav data...')
load(sprintf('Data_OOD/%s%d/%s_behavMeas.mat', subjName, nblocks, subjName))
fprintf('DONE\n')

load('Data_OOD/params_RC.mat')

% check if the iPRS of each pair are the same
for ipair = 1:max(dataMatrix(:, 8))
    ii = find(dataMatrix(:, 8)==ipair); % Col#8 is ipair
    iPrsA = dataMatrix(dataMatrix(:, 1)==ii(1), 6); % % Col#1 is itrial_allT, Col#6 is iPRS
    iPrsB = dataMatrix(dataMatrix(:, 1)==ii(2), 6);
    if iPrsA ~= iPrsA, fprintf('when ipair=%d, iPrsA=%d iPrsB=%d\n', ii, iPrsA, iPrsB), end
end

%% load source energy
nameEnergySource = sprintf('Data_OOD/%s%d/%s_energy_%d_%d.mat', subjName, nblocks, subjName, nORI, nSF);

fprintf('Loading the source energy ...')
load(nameEnergySource)
fprintf('DONE\n\n')

nAllTrials = size(e2D_allT, 1);
nAllTrials_perLoc = nAllTrials/nLoc5;

%% split and standardize energy
if iLocComb < 6, iLoc_all = iLocComb;
else, switch iLocComb , case 6, iLoc_all = [2,4]; case 7, iLoc_all = [3,5]; case 8, iLoc_all = 2:5; end
end

ntrain = round(nAllTrials_perLoc * ratio_train/2)*2; % to make sure ntrain is an even number
ntest = nAllTrials_perLoc - ntrain;

%% empty containers
data_allB = cell(ni, 1);
data_metrics_allB = nan(ni, nmetrics);
thresh_est_allB = nan(ni, 1);
params_est_allB = nan(ni, nparams_model);
pred_metrics_allB = nan(ni, nmetrics);
nLL_allB = nan(ni,1);
IC_allB = nan(ni, 3); % AIC_ML, AICc_ML and BIC
h_allB = nan(ni, length(namesNormalityTest), 2); %
nResample_allB = nan(ni, 1);
p_est_allB = nan(ni, nparams2);
metrics2_test_allB = nan(ni, nmetrics);
metrics2_pred_allB = nan(ni, nmetrics);

fprintf('Running...\n')
for ii = 1:ni
    nResample = 0;
    flagResample = 1;
    while flagResample && nResample < minResample
        %%%%%%%%%%%
        fxn_resampleTrials
        %%%%%%%%%%%
        
        % standardize the energy for training set
        e2D_train_norm = normEnergy(e2D_train, cst_train, iPRS_train);
        
        % obtain the kernels from the trained set of energy
        [kernel2D] = SX_sim07_RC(filtersSF_all, filtersOri_all, e2D_train_norm, resp_train);
        
        % derive the metrics for the test set
        pHit = mean(resp_test & iPRS_test)*2;
        pFA = mean(resp_test & (1-iPRS_test))*2;
        pC = (pHit+1-pFA)/2;
        [d,c] = SX_sim06_SDT(pHit, pFA);
        metrics_test = [d,c, [pC, pHit, pFA], nanmean(respC_test)]; % dprime, criterion, pC, pHit, pFA, pA, pA_PRS, pA_ABS
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % just estimate the thresh & mu/sigma given the data
        fxn_estMUSIGMA_fromMetrics
        p_est_allB(ii, :) = p_est;
        metrics2_test_allB(ii, :) = metrics_test;
        metrics2_pred_allB(ii, :) = metrics_pred;
        %%%%%%%%%%%%%%%%%%%%%%%%
        
        %  calculate the internal variable
        clear data % as the size of IV differs across subjects
        
        % template
        [template] = fxn_getTemplate(kernel2D, templateType);
        
        if iModelA == 3
            e2D_test = randn(ntest*length(iLoc_all), nORI, nSF) * std(e2D_test(:)) + mean(e2D_test(:));
        end
        e2D_test = e2D_test - min(e2D_test(:)) + eps;
        [IV_PRS_raw, IV_ABS_raw, maxPRS_allT, maxABS_allT] = fxn_getIV_v2(iModelA, e2D_test, iPRS_test, ...
            convolveType, IVType, template, ORI_bound);
        nPRS = length(IV_PRS_raw);
        nABS = length(IV_ABS_raw);
        
        if minResample == 1
            IV_PRS = IV_PRS_raw;
            IV_ABS = IV_ABS_raw;
        else
            if min(IV_PRS_raw)<=0, error('There is neg. values in IV PRS'); end
            if min(IV_ABS_raw)<=0, error('There is neg. values in IV ABS'); end
            IV_trans = boxcox([IV_PRS_raw; IV_ABS_raw]);
            IV_PRS = IV_trans(1:nPRS);
            IV_ABS = IV_trans(nPRS+1: end); assert(length(IV_ABS) == nABS)
%             [IV_PRS, lambda_PRS] = boxcox(IV_PRS_raw);
%             [IV_ABS, lambda_ABS] = boxcox(IV_ABS_raw);
        end
        
        % plot IV and the iMax
        plotFlag = 0; quickPlot_IV
        
        flagResample = swtest(IV_PRS) | swtest(IV_ABS);
        nResample = nResample+1;
        fprintf('*')
    end % while statement
    nResample_allB(ii)= nResample;
    
    %% compile data
    ndata = length(IV_PRS) + length(IV_ABS);
    data.ndata = ndata;
    data.IV{1} = IV_PRS; % before adding internal noise
    data.IV{2} = IV_ABS;
    data.imax{1} = maxPRS_allT;
    data.imax{2} = maxABS_allT;
    data.resp{1} = resp_test(boolean(iPRS_test));
    data.resp{2} = resp_test(boolean(1-iPRS_test));
    data.metrics_sim = metrics_test;     % keep the field name 'metrics_sim'!!
    data.cst = cst_test;
    data.iPRS = iPRS_test;
    data.respC = respC_test;
    data_allB{ii} = data;
    data_metrics_allB(ii, :) = metrics_test;
    
    % constrain the thresh
    % estimate the thresh given the measured criterion, to constrain the range of thresh when estimating params
    data.noisy_mu_PRS = median(IV_PRS);
    data.noisy_mu_ABS = median(IV_ABS);
    data.noisy_sigma_PRS = std(IV_PRS);
    data.noisy_sigma_ABS = std(IV_ABS);
    fxn_constrainThresh
    thresh_est_allB(ii) = thresh_est;
    
    %%%%%%%
    %   Fitting  %
    %%%%%%%
    fxn_estParams = @(params) fxn_getError_v2(iModelB, flagIncludePA, params, data);
    problem_ML = createOptimProblem('fmincon','objective', fxn_estParams,'x0', params0, 'lb', params_lb, 'ub', params_ub, 'options',options);
    ms_ML = MultiStart('StartPointsToRun', 'bounds', 'Display', 'off', 'UseParallel', 1);
    [params_est, nLL] = run(ms_ML, problem_ML, nrep);
    params_est_allB(ii, :) = params_est;
    nLL_allB(ii) = nLL;
    IC_allB(ii, :) = [getAIC_ML(nLL, nparams_model, ndata), ...
        getAICc_ML(nLL, nparams_model, ndata), ...
        getBIC(nLL, nparams_model, ndata)];
    
    %%%%%%%%%
    %   Prediction  %
    %%%%%%%%%
    pred = PR_pred(iModelB, flagIncludePA, params_est, data);
    pred_metrics_allB(ii, :) = pred.metrics;
    
    if find(ii == indCountdown1), fprintf('%d\n', indCountdown2(round(ii/(ni/nCountdown )))), end
    
end % end of ii

fprintf('\n\nALL iterations DONE\n')

%%%%%%%%%
switch iLocComb, case 1, color_ = 'r'; case 2, color_ = 'b'; end
% ModelPlot_local_perSubj
%%%%%%%%%

save(nameFileModelIDVD, '*_allB', 'names*', 'flag*', 'i*', 'minResample', 'ratio_train', 'ORI_bound', '*Type')
