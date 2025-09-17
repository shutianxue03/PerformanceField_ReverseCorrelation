function [k, i, margORI, margSF, R2, R2_Tjur, p, pCat, sep] = SX_RC6_kernel_parfor(e3D, dataMatrix, filtersSF_all, filtersOri_all)
% SX_RC6_kernel_parfor: Computes regression kernels and statistical metrics
%
% This function estimates regression kernels and various statistical metrics
% (e.g., R-squared, p-values) across different conditions (PRS, ABS, and BOTH)
% for orientation and spatial frequency channels. 

% INPUTS:
%   e3D              - 3D array of energy values (trials x orientations x SFs)
%   dataMatrix       - Matrix containing trial data, where column 6 indicates PRS/ABS
%   filtersSF_all    - Filter settings for spatial frequencies (SFs)
%   filtersOri_all   - Filter settings for orientations (ORI)
%
% OUTPUTS:
%   k                - 2D regression kernel values for each type (nTypes x ORI x SF)
%   i                - Intercept values for each regression (nTypes x ORI x SF)
%   margORI          - Marginalized kernel values across orientations (nTypes x ORI)
%   margSF           - Marginalized kernel values across spatial frequencies (nTypes x SF)
%   R2               - R-squared values for model fits (nTypes x ORI x SF)
%   R2_Tjur          - Tjur R-squared values for model fit quality (nTypes x ORI x SF)
%   p                - p-values for intercept and slope (nTypes x ORI x SF x 2)
%   pCat             - Categorical p-values for each type (nTypes x ORI x SF)
%   sep              - Separability check metric for each type (nTypes x 1)

% Define the number of trial types (1=PRS, 2=ABS, 3=BOTH)
nTypes = 3;
[ntrials, nORI, nSF] = size(e3D);

%% Initialize empty containers for outputs
k = nan(nTypes, nORI, nSF);              % Kernel values for each type
i = k;                                   % Intercepts for each type
R2 = k;                                  % R-squared values for each type
R2_Tjur = k;                             % Tjur R-squared values
margORI = nan(nTypes, nORI);             % Marginal values over ORI
margSF = nan(nTypes, nSF);               % Marginal values over SF
p = nan(nTypes, nORI, nSF, 2);           % p-values for intercept and slope
pCat = k;                                % Categorical p-values
sep = nan(nTypes, 1);                    % Separability metric

% Loop through each trial type (PRS, ABS, BOTH)
for iType = 1:nTypes
    % Define trial index based on type
    switch iType
        case 1, iPRS = dataMatrix(:, 6) == 1;  % PRS trials
        case 2, iPRS = dataMatrix(:, 6) == 0;  % ABS trials
        case 3, iPRS = 1:ntrials;              % BOTH (all trials)
    end
    
    % Perform regression and compute metrics using SX_sim07_RC function
    [kernel2D, intercept2D, R2_2D, R2_Tjur_, pValues_, pCat_] = SX_sim07_RC(filtersSF_all, filtersOri_all, e3D(iPRS, :, :), dataMatrix(iPRS, 9));
    
    % Store results in corresponding output containers
    k(iType, :, :) = kernel2D;
    i(iType, :, :) = intercept2D;
    R2(iType, :, :) = R2_2D;
    R2_Tjur(iType, :, :) = R2_Tjur_;
    p(iType, :, :, :) = pValues_;
    pCat(iType, :, :) = pCat_;
    
    % Compute separability metric for the current type
    sep(iType) = getSeparability(kernel2D);
    
    % Calculate marginalization across orientations and spatial frequencies
    margORI(iType, :) = squeeze(mean(kernel2D, 2));  % Average across SF
    margSF(iType, :) = squeeze(mean(kernel2D, 1));   % Average across ORI
end % End iType loop
