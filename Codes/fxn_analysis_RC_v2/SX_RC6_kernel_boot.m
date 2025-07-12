tic
%% boostrapping
% energy and resp of each type at each loc are randomly sampled (with the
% same index) WITH repetition, the fitted by logistic regression to get teh
% kernel

CI_level = .68; % confidence interval (plotted as errorbar)  (Fernandez, 2021)

if nB>1, fprintf(' boostrapping on kernel started (nB = %d):\n', nB), end

%% create empty containers
% typeInd = 2; exp_createContainers
kernels2D_allB = nan(nB, ntypes, nLoc, nfiltersOri, nfiltersSF); % 2D data
intercept_allB = kernels2D_allB;
R2_allB = nan(nB, ntypes, nLoc, nfiltersOri, nfiltersSF);
R2_Tjur_allB = R2_allB;
pValues_allB = nan(nB, ntypes, nLoc, nfiltersOri, nfiltersSF, 2); % 2 = p of intercept+p of slope
seprblity_kernel_allB = nan(nB, ntypes, nLoc); % p values of check separability
nAllTrials = ntrials*2;

%%
for iLoc = 1:nLoc
    if nB > 1, fprintf('    - %s:', namesLoc2D{iLoc}), end
    e_iLoc = energy_norm_perLoc{iLoc};
    resp_iLoc =  squeeze(data_both(:, iLoc, 2));
    
    for itype = 1:3   % PRS, ABS, BOTH
        
        
        
        energy = squeeze(e_iLoc(trialInd, :, :)); % resample energy WITH repetition
        resp = resp_iLoc(trialInd); % resample resp WITH repetition
        [kernel2D, yfit, R2, R2_Tjur, pValues, intercept] = SX_sim07_RC(filtersSF_all, filtersOri_all, energy, resp);
        
        % 1. save data
        intercept_allB(iB, itype, iLoc, :, :) = intercept;
        kernels2D_allB(iB, itype, iLoc, :, :) = kernel2D;
        R2_allB(iB, itype, iLoc, :, :) = R2;
        R2_Tjur_allB(iB, itype, iLoc, :, :) = R2_Tjur;
        pValues_allB(iB, itype, iLoc, :, :, :) = pValues;
        
        % 2. separability
        seprblity_kernel_allB(iB, itype, iLoc) = getSeparability(kernel2D);
    end
    if nB>1, fprintf(' %s ', namesType{itype}), end
end

%% get the average and CI
kernels2D = squeeze(median(kernels2D_allB,1));
intercept = squeeze(median(intercept_allB,1));
R2 = squeeze(median(R2_allB,1));
R2_Tjur = squeeze(median(R2_Tjur_allB,1));
pValues = squeeze(median(pValues_allB,1));
seprblity_kernel = squeeze(median(seprblity_kernel_allB,1));

kernels2D_lb = zeros(ntypes, nLoc, nfiltersOri, nfiltersSF);
kernels2D_ub = kernels2D_lb;
seprblity_kernel_lb = zeros(ntypes, nLoc);
seprblity_kernel_ub = seprblity_kernel_lb;

if nB > 1
    for itype = 1:3   % PRS, ABS, BOTH
        for iLoc = 1:nLoc
            kk = squeeze(kernels2D_allB(:, itype, iLoc, :, :));
            kernels2D_lb(itype, iLoc, :, :) = squeeze(quantile(kk,  .5-CI_level/2));
            kernels2D_ub(itype, iLoc, :, :)  = squeeze(quantile(kk,  .5+CI_level/2));
            
            ss = squeeze(seprblity_kernel_allB(:, itype, iLoc));
            seprblity_kernel_lb(itype, iLoc)  = squeeze(quantile(ss,  .5-CI_level/2));
            seprblity_kernel_ub(itype, iLoc)  = squeeze(quantile(ss,  .5+CI_level/2));
        end
    end
end
toc
%%
if energyFlag == 1 % only save var given using patch to avoid overwriting the saved files
    save(RCMatName, 'rHit_tgt', 'rFA_tgt', 'rHit_allFilters', 'rFA_allFilters', ...
        'ebin_prs_tgt', 'ebin_abs_tgt', 'ebin_prs_allFilters', 'ebin_abs_allFilters', ...
        'kernels2D', 'seprblity_kernel', 'R2', 'R2_Tjur', 'pValues', ...
        'kernels2D_lb', 'kernels2D_ub', 'seprblity_kernel_lb', 'seprblity_kernel_ub')
end

