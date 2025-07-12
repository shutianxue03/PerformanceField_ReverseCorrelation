
%%
itype = 3;
nrep = 20;
sigma_ext = .2;

%% extract pA & pC for each unik cst
nblocks_unik_perLoc = cell(nLoc5, 1);
cst_unik_perLoc = nblocks_unik_perLoc;
pA_unik_perLoc = nblocks_unik_perLoc;
pC_unik_perLoc = nblocks_unik_perLoc;
prop_perLoc = nblocks_unik_perLoc;

for iLoc = 1:nLoc5
    % do NOT % collapse all pA/pC of thesame cst
    nblocks = size(cst_perSess_perLoc, 1);
    
    % empty containers
    %     nblocks_perCST = nan(nblocks, 1); % number of blocks per contrast
    %     pA_perCST = nblocks_perCST;
    %     pC_perCST = nblocks_perCST;
    %
    %     % compile pA/C of each unik cst
    %     for ibk = 1:nblocks
    %
    %         ind_unik = cst_unik == cst_perSess_perLoc(:, iLoc);
    %         nblocks_perCST(ibk) = sum(ind_unik); % number of blocks per contrast
    %         pA_perCST(ibk) = mean(pA_perSess_perLoc(itype, ind_unik, iLoc));
    %         pC_perCST(ibk) = mean(pC_perSess_perLoc(itype, ind_unik, iLoc));
    %     end
    cst_unik_perLoc{iLoc} = cst_perSess_perLoc(:, iLoc);
    pA_unik_perLoc{iLoc} = squeeze(pA_perSess_perLoc(itype, :, iLoc))';
    pC_unik_perLoc{iLoc} = squeeze(pC_perSess_perLoc(itype, :, iLoc))';
    prop_perLoc{iLoc} = ones(nblocks, 1);
    
    % collapse all pA/pC of the same cst
    [cst_unik_all, ~, ~] = unique(cst_perSess_perLoc(:, iLoc)); % unique contrasts used at this loc
    ncst = length(cst_unik_all); % number of unique contrasts used at this loc
    %
    %     % empty containers
    nblocks_perCST = nan(ncst, 1); % number of blocks per contrast
    %     pA_perCST = nblocks_perCST;
    %     pC_perCST = nblocks_perCST;
    %
    %     % compile pA/C of each unik cst
    for icst_unik = 1:ncst
        cst_unik = cst_unik_all(icst_unik);
        ind_unik = cst_unik == cst_perSess_perLoc(:, iLoc);
        nblocks_perCST(icst_unik) = sum(ind_unik); % number of blocks per contrast
        %         pA_perCST(icst_unik) = mean(pA_perSess_perLoc(itype, ind_unik, iLoc));
        %         pC_perCST(icst_unik) = mean(pC_perSess_perLoc(itype, ind_unik, iLoc));
    end
    nblocks_unik_perLoc{iLoc} = nblocks_perCST; % unik stands for unique
    %     cst_unik_perLoc{iLoc} = cst_unik_all;
    %     pA_unik_perLoc{iLoc} = pA_perCST;
    %     pC_unik_perLoc{iLoc} = pC_perCST;
    %         prop_perLoc{iLoc} = nblocks_perCST/sum(nblocks_perCST);
end % end of iLoc

%% fitting
% Luzardo&Yeshurun, 2021
% pC vs. pA, estimate pC and pA by PTM
% range comes from Lu & Dosher 2008, Table 4 on p22
% the exact value comes from luzardo&yeshurun(2021)'s r code
% PTM
params0 = [2, 1, 1, .001, .5]; % (1) beta1, (2) beta2, (3) gamma, (4) sigma_add, (5) N_mul
params_lb =[0, 0, 2, 0, 0];
params_ub =[3, 3, 2.5, 3, 3];
namesParams = {'Beta1', 'Beta2', 'Gamma', 'Additive internal noise (sigma)','Multiplicative internal noise'};

% LAM
params0 = [2, 001, .5]; % (1) beta (2) sigma_add, (3) N_mul
params_lb =[0,  0, 0];
params_ub =[3, 3, 3];
namesParams = {'Beta', 'Additive internal noise (sigma)','Multiplicative internal noise'};

nparams = length(params_ub);
options = optimoptions('fmincon','MaxIterations',5000,'Display','off');
fileName_PTM = sprintf('Data/%s_PTM.mat', subjName);
fileDir_PTM = dir(fileName_PTM);

if isempty(fileDir_PTM)
    params_est_unik_perLoc = nan(nLoc5, nparams);
    pA_pred_unik_perLoc = cell(nLoc5, 1);
    pC_pred_unik_perLoc = pA_pred_unik_perLoc;
    
    clc
    fprintf('Fitting started...\n')
    for iLoc = 1:nLoc5
        %  extract data
        cst_unik_all = cst_unik_perLoc{iLoc};
        pA_perCST = pA_unik_perLoc{iLoc};
        pC_perCST = pC_unik_perLoc{iLoc};
        prop = prop_perLoc{iLoc};
        ncst = length(cst_unik_all);
        
        % fit and predict (unik cst)
        [params_est_unik, pA_pred_unik, pC_pred_unik] = fxn_FitPred(pA_perCST, pC_perCST, cst_unik_all, sigma_ext, 1, params0, params_lb, params_ub, options, nrep);
        params_est_unik_perLoc(iLoc, :) = params_est_unik;
        pA_pred_unik_perLoc{iLoc} = pA_pred_unik;
        pC_pred_unik_perLoc{iLoc} = pC_pred_unik;
        
        fprintf(sprintf('Loc%d done\n', iLoc))
    end % end of iLoc
    save(fileName_PTM, 'params_est_unik_perLoc', 'pA_pred_unik_perLoc', 'pC_pred_unik_perLoc')
else
    load(fileName_PTM)
end


%% plot data and prediction
ncomp = 3;
figure('Position', [0 0 1500 ncomp*300])

for iLoc = 1:nLoc5
    nblocks_perCST = nblocks_unik_perLoc{iLoc};
    cst_unik_all = cst_unik_perLoc{iLoc};
    pA_perCST = pA_unik_perLoc{iLoc};
    pC_perCST = pC_unik_perLoc{iLoc};
    prop = prop_perLoc{iLoc};
    params_est_unik = params_est_unik_perLoc(iLoc, :);
    pA_pred = pA_pred_unik_perLoc{iLoc};
    pC_pred = pC_pred_unik_perLoc{iLoc};
    
    ncst = length(nblocks_perCST);
    for icomp = 1:ncomp
        subplot(ncomp, nLoc5, (icomp-1)*nLoc5+iLoc), hold on, grid on, box on
        switch icomp
            case 1,  x = pA_perCST; y = pC_perCST; x_pred = pA_pred; y_pred = pC_pred; xlabel_ = 'pA'; ylabel_ = 'pC';
            case 2, x = cst_unik_all; y = pC_perCST; x_pred = x;            y_pred = pC_pred; xlabel_ = 'signal cst'; ylabel_ = 'pC';
            case 3, x = cst_unik_all; y = pA_perCST; x_pred = x;            y_pred = pA_pred; xlabel_ = 'signal cst'; ylabel_ = 'pA';
        end
        %%%%%%%%
        % plot DATA %
        %%%%%%%%
        for ibk = 1:ncst
            % plot data points AVERAGED across CST levels
            markerSize = prop(ibk)*5+5;
            plot(x(ibk), y(ibk), 'o', 'color', colors_comb(iLoc, :), 'markersize', markerSize)
        end
        
        %%%%%%%%%%%
        %  plot predictions  %
        %%%%%%%%%%%
        [x_pred_sorted, i] = sort(x_pred);
        plot(x_pred_sorted, y_pred(i), 'k-', 'linewidth', 2)
        
        if icomp == 1
            xlim([.5, 1]), ylim([.5, 1])
            title(namesLoc2D{iLoc})
        end
        xlabel(xlabel_)
        ylabel(ylabel_)
    end % end of icomp
end % end of iLoc
sgtitle('Predicted by PTM')

%% plot estimated params vs. pA
figure('Position', [0 0 1200 200])
for iparam = 1:nparams
    subplot(1,nparams,iparam), hold on
    for iLoc = 1:nLoc5
        pA_perCST = pA_unik_perLoc{iLoc};
        params_est_unik = params_est_unik_perLoc(iLoc, :);
        errorbar(params_est_unik(iparam), mean(pA_perCST), std(pA_perCST), 'o', 'color', colors5Loc(iLoc, :), 'Capsize', 0)
    end
    xlabel(namesParams{iparam})
    ylabel('pA')
    ylim([.5, 1])
    
end
set(findall(gcf, '-property', 'fontsize'), 'fontsize',12)
sgtitle('Estimated parameters')


%%
function [params_est, pA_pred_perC, pC_pred_perC] = fxn_FitPred(pA, pC, cst, sigma_ext, prop, params0, params_lb, params_ub, options, nrep)
ncst = length(cst);

% fit data (unik cst)
fxn_estParams = @(params) fxn_getError_IN(params, pA, pC, cst, sigma_ext, prop);
problem_ML = createOptimProblem('fmincon','objective', fxn_estParams,'x0', params0, 'lb', params_lb, 'ub', params_ub, 'options',options);
ms_ML = MultiStart('StartPointsToRun', 'bounds', 'Display', 'off', 'UseParallel', 1);
[params_est, fval, exitflag, output, solutions] = run(ms_ML, problem_ML, nrep);
%     quickPlot_globalMin(solutions, nparams, params_est, params_lb,params_ub, fval, nrep)

% make predictions
pA_pred_perC = nan(ncst,1);
pC_pred_perC = pA_pred_perC;
for icst = 1:ncst
    %     [pA_pred, pC_pred] = fxn_pred_pApC_PTM(params_est, cst(icst), sigma_ext);
    [pA_pred, pC_pred] = fxn_pred_pApC_LAM(params_est, cst(icst), sigma_ext);
    pA_pred_perC(icst) = pA_pred;
    pC_pred_perC(icst) = pC_pred;
end
end


%% helper funcion - get error
function error = fxn_getError_IN(params, pA_data, pC_data, cst, sigma_ext, prop)

ncst = length(pA_data);
pA_pred_all = nan(ncst, 1);
pC_pred_all = pA_pred_all;
for icst = 1:ncst
    %     [pC_pred, pA_pred] = fxn_pred_pApC_PTM(params, cst(icst), sigma_ext);
    [pC_pred, pA_pred] = fxn_pred_pApC_LAM(params, cst(icst), sigma_ext);
    pA_pred_all(icst)= pC_pred;
    pC_pred_all(icst)= pA_pred;
end

% error = sumsqr((pA_pred_all-pA_data).*prop) + sumsqr((pC_pred_all-pC_data).*prop);
error = sumsqr(pA_pred_all-pA_data) + sumsqr(pC_pred_all-pC_data);
end

%% helper funcion - predict pA and pC
function [pC,pA] = fxn_pred_pApC_PTM(params, c, sigma_ext)

beta1 = params(1);         % gain factor
beta2 = params(2);        % gain factor
gamma = params(3);      % power factor
sigma_add = params(4);  % additive internal noise
N_mul = params(5);        % multiplicative internal noise

pC_PDF_sigma = sqrt(sigma_ext^(2*gamma)+N_mul^2*(sigma_ext^(2*gamma)+(beta2*c)^(2*gamma))+sigma_add^2);
pC_CDF_sigma = sqrt(sigma_ext^(2*gamma)+N_mul^2*sigma_ext^(2*gamma)+sigma_add^2);
pA_CDF_sigma = sqrt(N_mul^2*(2*sigma_ext^(2*gamma)+(beta2*c)^(2*gamma)) + 2*sigma_add^2);

pC_PDF = @(x) normpdf(x - (beta1*c)^gamma, 0, pC_PDF_sigma);
pC_CDF = @(x) normcdf(x, 0, pC_CDF_sigma);
pA_PDF = @(x) normpdf(x-(beta1*c)^gamma, 0, sqrt(2)*sigma_ext^gamma);
pA_CDF = @(x) normcdf(x, 0, pA_CDF_sigma);

fxn_pC = @(x) pC_PDF(x).*pC_CDF(x);
fxn_pA = @(x) pA_PDF(x) .* (pA_CDF(x) .^ 2 + (1-pA_CDF(x)) .^ 2);
pC = integral(fxn_pC, -inf, inf);
pA = integral(fxn_pA, -inf, inf);
end

%%
function [pC,pA] = fxn_pred_pApC_LAM(params, c, sigma_ext)

beta = params(1);         % gain factor
sigma_add = params(2);  % additive internal noise
N_mul = params(3);        % multiplicative internal noise

pA_CDF_sigma = sqrt(2*(sigma_add^2+N_mul^2*sigma_ext^2));
pC_PDF = @(x) normpdf(x - (beta*c)/sqrt((1+N_mul^2)*sigma_ext^2 + sigma_add^2), 0,1);
pC_CDF = @(x) normcdf(x, 0, 1);
pA_PDF = @(x) normpdf(x-beta*c, 0, sqrt(2)*sigma_ext);
pA_CDF = @(x) normcdf(x, 0, pA_CDF_sigma);

fxn_pC = @(x) pC_PDF(x).*pC_CDF(x);
fxn_pA = @(x) pA_PDF(x) .* (pA_CDF(x) .^ 2 + (1-pA_CDF(x)) .^ 2);
pC = integral(fxn_pC, -inf, inf);
pA = integral(fxn_pA, -inf, inf);
end


%% helper fxn - quickPlot
% plot all local solution at eahc iteration
function quickPlot_globalMin(solutions, nparams, params_est, params_lb,params_ub, fval, nrep)
nS = length(solutions);
params_est_allS = reshape([solutions.X],  nparams, nS); % confirmed 100% correct
params0_allS = nan(nparams, nS); for is = 1:nS, params0_allS(:, is) = solutions(is).X0{1}; end
paramsFVAL_allS = [solutions.Fval];

buffer = .1;
figure('Position', [3000 0 1000 1000])
for iparam = 1:nparams
    subplot(nparams+1, 1, iparam), hold on
    plot((1:nS)-buffer, params_est_allS(iparam, :), 'o')
    plot((1:nS)+buffer, params0_allS(iparam, :), '+')
    yline(params_est(iparam), 'k-');
    if iparam == 1, legend({'estimation', 'start'}, 'Location', 'best'), end
    xlim([0, nS+1])
    ylim([params_lb(iparam), params_ub(iparam)])
end
subplot(nparams+1, 1, nparams+1), hold on
plot(1:nS, paramsFVAL_allS, 'ko')
yline(fval, 'k-');
xlim([0, nS+1])

sgtitle(sprintf('%d/%d reps worked out', nS, nrep))

figure('Position', [4000 0 100 1000])
for iparam = 1:nparams
    subplot(nparams+1, 1, iparam), hold on
    histogram(params_est_allS(iparam, :), 5, 'DisplayStyle', 'Stair')
    xline(params_est(iparam), 'r-');
    xlim([params_lb(iparam), params_ub(iparam)])
end
end
