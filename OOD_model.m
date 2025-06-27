function OOD_model(isubj, iLocComb, iModelA, iModelB, convolveType, IVType, ni, flag_PatchMode, flag_standEnergy)

%% INPUT
%      isubj: index of  subj
%      iLocComb: 1=Fovea, 8=periF(6 deg ecc), 6=HM, 7=VM, 5=LVM, 3=UVM
%      iModelA: 1=core model, 2=randomize template, 3=randomize energy,
%                     4=randomize both template and energy; 5=use IO template
%      iModelB: 4=do NOT estimate internal noise, 1=estimate both induced & constant noise, 2=estimate only constant noise, 3=estimate only induced noise
%      convolveType: 1=dot multiply template and energy, 2=convolve
%      IVType: 1=sum, 2=max
%      ni: number of iterations
%      flag_PatchMode: when estimating templates, 1=energy derives from target; 2=noise patch
%      flag_standEnergy: 1=standardize energy, 0=do NOT
%      flagConstrainThresh: 0= not constrain threshold

clc
close all
time_start = datetime('now')
warning off
format compact

addpath(genpath('Data_OOD'))
addpath(genpath('Data_model'))
addpath(genpath('fxn_model'))
addpath(genpath('fxn_RCplot'))
addpath(genpath('fxn_analysis_RC_v2'))
load('Data_OOD/params.mat')

flag_cutMapping=0;
SX_RC1_setting

%% general params
nrep = 15; % number of repetitions to look for global minimum in fmincon
templateType = 1; % (1) mirrored kernel (2) positive kernel (3) reconstructed kernel
itype_template = 2; % 1=estimate template from PRS trials, ABS trials, or BOTH trials
iSess_start = 6; % from which session data is taken into account
minResample = 1;
flagConstrainThresh=0;

switch iModelB
    case 4, params0=[nan];params_lb=[nan];params_ub=[nan];
    case 1, params0 = [1, .5, nan]; params_lb = [1e-5, 1e-5, nan]; params_ub = [5, 2, nan];
    case 2, params0 = [.5, nan]; params_lb = [1e-5, nan]; params_ub = [5, nan];
    case 3, params0 = [1, nan]; params_lb = [1e-5, nan]; params_ub = [5, nan];
end

%%
SX_RC1_setting

flag_plotIVsDist = 0;
ratio_train = 3/4; % proportion of trials in the training set (to derive 2D kernels); others are in the test set (to estimate params)
ORI_bound = [10, 20]; % -40 and 40

getAIC_ML = @(nLL, nparamsAll, ndata) 2 * nLL + 2*nparamsAll;
% getAIC_LS = @(err, nparamsAll, ndata) ndata * log(err/ndata) + 2*nparamsAll;
getAICc_ML = @(nLL, nparamsAll, ndata) 2 * nLL + 2*nparamsAll + 2*nparamsAll*(nparamsAll+1)/(ndata-nparamsAll-1);
% getAICc_LS = @(err, nparamsAll, ndata) ndata * log(err/ndata) + 2*nparamsAll + 2*nparamsAll*(nparamsAll+1)/(ndata-nparamsAll-1);
getBIC = @(nLL, nparamsAll, ndata) 2*nLL + nparamsAll * log(ndata);

options = optimoptions('fmincon','MaxIterations', 1e4, 'Display','off');

markers_allSubj = {'o', 's', 'd', '^','v',  '<', '>','p', 'h', '+', 'x', '-'}; % for each subj

%% names
% model_setting
subjList =             {'YK', 'SP', 'SX', 'LS', 'RE', 'MD', 'AS', 'HL', 'FH', 'HA', 'CS', 'DT', 'DU', 'RC','SR'};
nblocks_allSubj = [200, 240, 210, 240, 220, 220, 215, 220, 205, 210, 205, 205, 205, 205, 195];

namesParamsModel = namesParamsModel_all{iModelB};
nparams_model = length(namesParamsModel); % number of params to fit; currently: alpha and thresh
assert(nparams_model == length(params0))

subjName = subjList{isubj};
nblocks = nblocks_allSubj(isubj);

%% file name & print
% nameFolder_model = sprintf('Data_NOM/ORI%dSF%d', nORI, nSF);
% if isempty(dir(nameFolder_model)), mkdir(nameFolder_model), end
% if flag_PatchMode == 1
%     nameFileModelIDVD = sprintf('%s/%s_n%d_Loc%d_T_A%dB%d_min%d_thresh%d_conv%d_IV%d', ...
%         nameFolder_model, subjName, ni, iLocComb, iModelA, iModelB, minResample, flagConstrainThresh, convolveType, IVType);
% else
%     nameFileModelIDVD = sprintf('%s/%s_n%d_Loc%d_N%d_A%dB%d_min%d_thresh%d_conv%d_IV%d', ...
%         nameFolder_model, subjName, ni, iLocComb, flag_standEnergy, iModelA, iModelB, minResample, flagConstrainThresh, convolveType, IVType);
% end
nameFolder_NOM = sprintf('Data_NOM/ORI%dSF%d/%s/L%d', nORI, nSF, subjName, iLocComb);
if isempty(dir(nameFolder_NOM)), mkdir(nameFolder_NOM), end

nameFileModelIDVD = sprintf('%s/n%d_A%dB%d', nameFolder_NOM, ni, iModelA, iModelB);

fprintf('\n%s nblocks = %d [ORI%d SF%d]\n - ni = %d\n - Loc: %s\n - MODEL: [A%d] %s & [B%d] %s\n - Template (1=PRS, 2=ABS, 3=BOTH): %d\n - Min. resampling = %d\n - constrain threshold based on criterion: %d\n - Energy source: %d (1=TARGET, 2=NOISE)\n - Convolve type: %d\n - IV type: %d\n', ...
    subjName, nblocks, nORI, nSF, ni, namesLocComb{iLocComb}, iModelA, namesModelA{iModelA}, iModelB, namesModelB{iModelB}, itype_template, minResample, ...
    flagConstrainThresh, flag_PatchMode, convolveType, IVType)

%% load behav data
fprintf('\nLoading the behav data...')
% load(sprintf('Data_OOD/%s%d/%s_behavMeas.mat', subjName, nblocks, subjName))
load(sprintf('Data_OOD/%s%d/%s_behavMeas.mat', subjName, nblocks, subjName))
fprintf('DONE\n')

SX_RC1_setting
nmetrics=8; % the 9th one is CS in SX_RC1_setting, but no need for model

%% check if the iPRS of each pair are the same
for ipair = 1:max(dataMatrix(:, 8))
    ii = find(dataMatrix(:, 8)==ipair); % Col#8 is ipair
    iPrsA = dataMatrix(dataMatrix(:, 1)==ii(1), 6); % % Col#1 is itrial_allT, Col#6 is iPRS
    iPrsB = dataMatrix(dataMatrix(:, 1)==ii(2), 6);
    if iPrsA ~= iPrsA, fprintf('when ipair=%d, iPrsA=%d iPrsB=%d\n', ii, iPrsA, iPrsB), end
end

%% load source energy for the training group to derive the template
if flag_PatchMode == 1
    nameSourceEnergy = sprintf('Data_OOD/%s%d/%s_energy_T_%d_%d.mat', subjName, nblocks, subjName, nORI, nSF);
else
    nameSourceEnergy = sprintf('Data_OOD/%s%d/%s_energy_N_%d_%d.mat', subjName, nblocks, subjName, nORI, nSF);
end
fprintf('Loading the source energy (training set) ...')
load(nameSourceEnergy)
fprintf('DONE\n')
if flag_PatchMode == 1, e3D_allT_train = e3D_target_allT; else, e3D_allT_train = e3D_noise_allT; end

%% load source energy for the test group (must be energy derived from targets)
nameSourceEnergy = sprintf('Data_OOD/%s%d/%s_energy_T_%d_%d.mat', subjName, nblocks, subjName, nORI, nSF);
fprintf('Loading the source energy (test set) ...')
load(nameSourceEnergy)
fprintf('DONE\n\n')
e3D_allT_test = e3D_target_allT;

%% trial number and loc index
nAllTrials = size(e3D_allT_train, 1); % number of total trials collected

nSess = nAllTrials/500; % nunber of sessions per single loc
iSess_select = iSess_start:nSess; % the index of selected sessions
nAllTrials = nAllTrials * length(iSess_select)/nSess;
ntrials_perSingleLoc = nAllTrials/nLoc5; % number of trial per single location, i.e., fovea/L/R/UVM/LVM
nLocRep = 4;
ntrials_perCompLoc = ntrials_perSingleLoc * nLocRep; % number of trial per location for comparison, i.e., fovea/P/HM/VM/LVM/UVM

% to equate the number of trials per loc, single location got
% doubled (HM vs. VM) even quandrupled (F vs. P)
if iLocComb < 6, iLoc_all = ones(1,nLocRep)*iLocComb;
else, switch iLocComb , case 6, iLoc_all = [2,2,4,4]; case 7, iLoc_all = [3,3,5,5]; case 8, iLoc_all = 2:5; end
end
assert(length(iLoc_all) == nLocRep)

% no. trials for training and testing group per location to be compared (F, P, HM, VM, L, R, UVM, LVM)
ntrain_perSingle = ntrials_perSingleLoc * ratio_train; if ntrain_perSingle<1e3, error('WARINING: training set size is < 1K, might damage RC quality\n'), end% since ratio is 3/4, ntrain is 100% certain an integer
ntest_perSingle = ntrials_perSingleLoc * (1-ratio_train);
if rem(ntrain_perSingle, 2), ntrain_perSingle = ntrain_perSingle+1;ntest_perSingle = ntest_perSingle-1; end
ntrain_perComp = ntrain_perSingle * nLocRep; % since ratio is 3/4, ntrain is 100% certain an integer
ntest_perComp = ntest_perSingle * nLocRep;

%% empty containers
data_allB = cell(ni, 1);
data_metrics_allB = nan(ni, nmetrics);
thresh_est_allB = nan(ni, 1);
params_est_allB = nan(ni, nparams_model);
pred_metrics_allB = nan(ni, nmetrics);
nLL_allB = nan(ni,1);
IC_allB = nan(ni, 3); % AIC_ML, AICc_ML and BIC
h_allB = nan(ni, length(namesNormalityTest), 2); % 2: PRS and ABS
nResample_allB = nan(ni, 1);
p_est_allB = nan(ni, nparams2);
metrics2_test_allB = nan(ni, nmetrics);
metrics2_pred_allB = nan(ni, nmetrics);

fprintf('Running...\n')
datetime('now')
for ii = 1:ni
%     nResample = 0;
%     flagResample = 1;
%     while flagResample && nResample < minResample
        
        %%%%%%%%%%%
        fxn_resampleTrials
        %%%%%%%%%%%
        e3D_test_true = e3D_test;
        
        % standardize the energy for the training set
        switch itype_template
            case 1
                e3D_train = e3D_train(iPRS_train, :, :);
                cst_train = cst_train(iPRS_train);
                resp_train = resp_train(boolean(iPRS_train));
                iPRS_train = iPRS_train(iPRS_train);
            case 2
                e3D_train = e3D_train(boolean(1-iPRS_train), :, :); % only use signal-ABS trials to derive the template
                cst_train = cst_train(boolean(1-iPRS_train));
                resp_train = resp_train(boolean(1-iPRS_train));
                iPRS_train = iPRS_train(boolean(1-iPRS_train));
        end
        
        if flag_standEnergy
            e3D_train_norm = normEnergy(e3D_train, cst_train, iPRS_train);
        else
            e3D_train_norm = e3D_train;
        end
        
        % obtain kernels from the training set
        
        kernel2D = SX_sim07_RC(filtersSF_all, filtersOri_all, e3D_train_norm, resp_train); % 5 seconds
%         figure, imagesc(kernel2D), axis square
        
        % derive the metrics for the test set
        pHit = mean(resp_test & iPRS_test)*2;
        pFA = mean(resp_test & (1-iPRS_test))*2;
        pC = (pHit+1-pFA)/2;
        [d,c] = SX_sim06_SDT(pHit, pFA);
        metrics_test = [d,c, [pC, pHit, pFA], nanmean(respC_test)]; % dprime, criterion, pC, pHit, pFA, pA, pA_PRS, pA_ABS
        metrics2_test_allB(ii, :) = metrics_test;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% directly estimate the thresh & mu/sigma from the data
        % (assuming normal distribution)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fxn_estMUSIGMA_fromMetrics
        p_est_allB(ii, :) = p_est; % thresh, mu_PRS, sigma_PRS, mu_ABS and sigma_ABS
        metrics2_pred_allB(ii, :) = metrics_pred;
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        %%  calculate internal variable
        %%%%%%%%%%%%%%%%%%%%%%%%
        clear data % as the size of IV differs across subjects
        
        %% 1. Template
        if iModelA== 5 % use IO template (the energy profile of the signal)
            if nORI==29, load('Data_OOD/signalEnergy_29_29.mat', 'signalEnergy'); end
            if nORI==37, load('Data_OOD/signalEnergy_37_37.mat', 'signalEnergy'); end
            template = signalEnergy;
        else
            [template] = fxn_getTemplate(kernel2D, templateType);
        end
        
        %% 2. Energy profile
        if (iModelA == 3) || (iModelA == 4) % randomize energy profile on each trial
            [ntrials, nORI, nSF] = size(e3D_test_true);
            e2D_test_inUse = nan(size(e3D_test_true));
            parfor itrial = 1:ntrials
                e2D_test_true = squeeze(e3D_test_true(itrial, :, :));
                rng('shuffle')
                indRand = randperm(nORI*nSF);
                e2D_test_true_v = e2D_test_true(:);
                e2D_test_v_rand = e2D_test_true_v(indRand);
                
                e2D_test_inUse_ = reshape(e2D_test_v_rand, nORI, nSF);
                assert(abs(sum(e2D_test_inUse_(:)) - sum(e2D_test_true(:))) < 1e-5)
                e2D_test_inUse(itrial, :, :) = e2D_test_inUse_;
                %                 figure, subplot(1,2,1), imagesc(e2D_test_true), axis square, subplot(1,2,2), imagesc(e2D_test_inUse_), axis square
            end % itrial
        else
            e2D_test_inUse = e3D_test_true;
        end
        
        %% 3. Internal variable (IV)
        % (the template is randomized, if needed)
        % around 13 seconds for SP
        [IV_PRS_raw, IV_ABS_raw, maxPRS_allT, maxABS_allT] = fxn_getIV_v2(iModelA, e2D_test_inUse, iPRS_test, ...
            convolveType, IVType, template, ORI_bound);
        nPRS = length(IV_PRS_raw);
        nABS = length(IV_ABS_raw);
        
        IV2 = [IV_PRS_raw; IV_ABS_raw]; IV2 = IV2 - min(IV2(:)) + eps;
        IV_trans = sqrt(IV2);
        IV_PRS = IV_trans(1:nPRS);
        IV_ABS = IV_trans(nPRS+1: end); assert(length(IV_ABS) == nABS)
        
        %% plot IV and the iMax
        IV2 = [IV_PRS_raw; IV_ABS_raw]; IV2 = IV2 - min(IV2(:)) + eps;
        %%%%%%%%%%%%%%%%%%%%%%
%         quickPlot_IV % output h_allB
        %%%%%%%%%%%%%%%%%%%%%%
        
%         flagResample = nan;%swtest(IV_PRS) | swtest(IV_ABS);
%         nResample = nResample+1;
%         fprintf('*')
%     end % while statement
%     nResample_allB(ii)= nResample;
   
    
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
    
    %%%%%%%%%%%%
    fxn_constrainThresh
    thresh_est_allB(ii) = thresh_est;
    %%%%%%%%%%%%
    
    %%%%%%%
    %   Fitting  %
    %%%%%%%
    fxn_estParams = @(params) fxn_getError_v4(iModelB, params, data);
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
    pred = PR_pred(iModelB, params_est, data);
    pred_metrics_allB(ii, :) = pred.metrics;
    
     % calculate human-computer consistency
        % resp % predicted trial-wise response
        % 
    if flag_plotIVsDist, quickPlot_noisyIV, end
    %     if find(ii == indCountdown1), fprintf('%d\n', indCountdown2(round(ii/(ni/nCountdown )))), end
    
    fprintf('%d/%d\n', ii, ni)
    % iModelA=1 (raw template & energy): <14 seconds
    % iModelA=2 (randomize template): <19 seconds
    % iModelA=3 (randomize energy): <21 seconds
%     nLL_4(iModelA) = median(nLL_allB);
    
end % end of ii

fprintf('\n\nALL iterations DONE\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_plotIVsDist, ModelPlot_local_perSubj, end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% do NOT save i* (in conflict when plotting)
save(nameFileModelIDVD, '*_allB', 'names*', 'flag*', 'minResample', 'ratio_train', 'ORI_bound', '*Type')
fprintf('Data saved\n')


time_end = datetime('now')

time_end - time_start
