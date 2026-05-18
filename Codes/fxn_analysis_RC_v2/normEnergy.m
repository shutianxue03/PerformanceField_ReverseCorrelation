function [e3D_norm, normStats] = normEnergy(e3D_allT, cst_allT, iPRS_allT)
% normEnergy: Standardizes energy values across signal conditions and contrast levels
%
% This function standardizes the energy values for each signal condition (e.g., PRS and ABS) 
% and contrast level, centering the data at 0 and scaling to a standard deviation (SD) of 1. 
% If the iPRS_allT array is not provided, the function assumes equal PRS and ABS trials.

% INPUTS:
%   e3D_allT    - 3D array of energy values for all trials
%   cst_allT    - Array indicating contrast levels for each trial
%   iPRS_allT   - Binary array indicating the signal condition (1 = PRS, 0 = ABS) [optional]

% OUTPUTS:
%   e3D_norm    - 3D array of normalized energy values
%   normStats   - Struct of normalization statistics with fields:
%                 .mu_byGroup       [nCst x 2 x nORI x nSF]
%                 .sigma_byGroup    [nCst x 2 x nORI x nSF]
%                 .mu_trial         [nTrials x nORI x nSF]
%                 .sigma_trial      [nTrials x nORI x nSF]
%                 .mu_global_3D     [1 x nORI x nSF]
%                 .sigma_global_2D  [nORI x nSF]
%                 .cst_unique       [1 x nCst]

% If only two arguments are given, assume equal PRS and ABS trials
if nargin == 2
    ntrials = size(e3D_allT, 1) / 2;
    iPRS_allT = [ones(ntrials, 1); zeros(ntrials, 1)];
    assert(length(iPRS_allT) == length(cst_allT), 'Mismatch in lengths of iPRS_allT and cst_allT.');
end

% Ensure logical PRS labels and column vectors for indexing.
iPRS_allT = logical(iPRS_allT(:));
cst_allT = cst_allT(:);

%% Preprocess energy values by contrast levels
e3D_norm = nan(size(e3D_allT));
[nTrials, nORI, nSF] = size(e3D_allT);

% Identify unique contrast levels in the data
cst_unik = unique(cst_allT);
cst_unik = cst_unik(:).';  % Ensure row vector for iteration
nCst = numel(cst_unik);

mu_byGroup = nan(nCst, 2, nORI, nSF);      % class index: 1=ABS, 2=PRS
sigma_byGroup = nan(nCst, 2, nORI, nSF);   % class index: 1=ABS, 2=PRS
mu_trial = nan(size(e3D_allT));
sigma_trial = nan(size(e3D_allT));
eps_sigma = 1e-8;

% Loop through each unique contrast level
for iCst = 1:nCst
    icst_unik = cst_unik(iCst);
    % Find indices for PRS and ABS trials at the current contrast level
    indPRS = (cst_allT == icst_unik & iPRS_allT);
    indABS = (cst_allT == icst_unik & ~iPRS_allT);

    if any(indABS)
        muABS = mean(e3D_allT(indABS, :, :), 1);
        sigmaABS = std(e3D_allT(indABS, :, :), [], 1);
        sigmaABS(~isfinite(sigmaABS) | sigmaABS < eps_sigma) = eps_sigma;
        e3D_norm(indABS, :, :) = bsxfun(@rdivide, bsxfun(@minus, e3D_allT(indABS, :, :), muABS), sigmaABS);
        mu_trial(indABS, :, :) = repmat(muABS, [sum(indABS), 1, 1]);
        sigma_trial(indABS, :, :) = repmat(sigmaABS, [sum(indABS), 1, 1]);
        mu_byGroup(iCst, 1, :, :) = muABS;
        sigma_byGroup(iCst, 1, :, :) = sigmaABS;
    end

    if any(indPRS)
        muPRS = mean(e3D_allT(indPRS, :, :), 1);
        sigmaPRS = std(e3D_allT(indPRS, :, :), [], 1);
        sigmaPRS(~isfinite(sigmaPRS) | sigmaPRS < eps_sigma) = eps_sigma;
        e3D_norm(indPRS, :, :) = bsxfun(@rdivide, bsxfun(@minus, e3D_allT(indPRS, :, :), muPRS), sigmaPRS);
        mu_trial(indPRS, :, :) = repmat(muPRS, [sum(indPRS), 1, 1]);
        sigma_trial(indPRS, :, :) = repmat(sigmaPRS, [sum(indPRS), 1, 1]);
        mu_byGroup(iCst, 2, :, :) = muPRS;
        sigma_byGroup(iCst, 2, :, :) = sigmaPRS;
    end
end

% Global summaries for compatibility with code paths expecting pooled stats.
X = reshape(e3D_allT, [nTrials, nORI * nSF]);
mu_global = mean(X, 1);
sigma_global = std(X, [], 1);
sigma_global(~isfinite(sigma_global) | sigma_global < eps_sigma) = eps_sigma;

normStats = struct();
normStats.mu_byGroup = mu_byGroup;
normStats.sigma_byGroup = sigma_byGroup;
normStats.mu_trial = mu_trial;
normStats.sigma_trial = sigma_trial;
normStats.mu_global_3D = reshape(mu_global, [1, nORI, nSF]);
normStats.sigma_global_2D = reshape(sigma_global, [nORI, nSF]);
normStats.cst_unique = cst_unik;

% 
% function e3D_norm = normEnergy(e3D_allT, cst_allT, iPRS_allT)
% 
% % standardize energy for each signal condition (PRS/ABS) and at each contrast 
% % to center at 0 and adjust SD to 1
% 
% if nargin == 2
%     ntrials = size(e3D_allT,1)/2;
%     iPRS_allT = [ones(ntrials,1); zeros(ntrials,1)];
%     assert(length(iPRS_allT) == length(cst_allT))
% end
% 
% %% preprocess based on cst levels (WITH itgt)
% e3D_norm = nan(size(e3D_allT));
% 
% % get unique cst
% cst_unik = unique(cst_allT);
% cst_unik = cst_unik(:).'; %reshape(cst_unique, 1, length(cst_unique));
% 
% % loop
% for icst_unik = cst_unik
%     indPRS = (cst_allT == icst_unik & iPRS_allT);
%     indABS = (cst_allT == icst_unik & ~iPRS_allT);
%     %     assert(sum(indPRS) == sum(indABS))
%     
%     % center at 0
%     e3D_allT(indPRS, :, :) = e3D_allT(indPRS, :, :) - mean(e3D_allT(indPRS, :, :), 1);
%     e3D_allT(indABS, :, :) = e3D_allT(indABS, :, :) - mean(e3D_allT(indABS, :, :), 1);
%     
%     % adjust SD to 1
%     e3D_norm(indPRS, :,:) = e3D_allT(indPRS, :, :) ./ std(e3D_allT(indPRS, :, :), [], 1);
%     e3D_norm(indABS, :,:) = e3D_allT(indABS, :, :) ./ std(e3D_allT(indABS, :, :), [], 1);
% %     energy_norm(indPRS, :,:) = energy(indPRS, :, :);
% %     energy_norm(indABS, :,:) = energy(indABS, :, :);
% end
% 
% %% AF's code
% % contra = unique(trialsMat(:,9)); %contrast
% %
% % for ii = 1:length(contra)
% %     ind1 = targettype==1 & trialsMat(:,9)==contra(ii);
% %     ind0 = targettype== 0 & trialsMat(:,9)==contra(ii);
% %     %subtract mean
% %     dX(:,:,ind1) = dX(:,:,ind1)-mean(dX(:,:,ind1),3); % the 3rd dim is itrial
% %     dX(:,:,ind0) = dX(:,:,ind0)-mean(dX(:,:, ind0),3);
% %     %divide by std
% %     dX(:,:,ind1) = dX(:,:,ind1)./std(dX(:,:,ind1),[],3);
% %     dX(:,:,ind0) = dX(:,:,ind0)./std(dX(:,:,ind0),[],3);
% % end
% 
% 
