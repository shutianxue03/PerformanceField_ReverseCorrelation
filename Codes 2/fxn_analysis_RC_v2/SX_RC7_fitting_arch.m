
function [data2D_perComb, marg_perComb, margPred_perComb, margParams_perComb, margR2_perComb] = SX_RC7_fitting_arch(data2D, nLoc, ntypes, nfiltersOri, nfiltersSF, namesLocComb, xaxis_all, namesModelParams, fitMode, ub_all, lb_all)

% make archive on Apr 21, 2022
% fit each function of the pair (e.g., Fovea vs. peri) individually
% while in the updated file, two functions are fitted together

mirrorFlag = 1; % mirror the 2D matrix over ORI=0
ncomb = 8;

imodelKernel_perF = [8,3]; 

%% 1. mirror the 2D image
if mirrorFlag == 1
    data2D_mirror = getMirror(data2D);
else
    data2D_mirror = data2D;
end

%% 3. combine locations
% look at my notes in notability
data2D_perComb = nan(ntypes, 8, nfiltersOri, nfiltersSF);
margOri_perComb = cell(length(namesLocComb), 1);
margSF_perComb = cell(length(namesLocComb), 1);

margDimOri = ndims(data2D_mirror);
margDimSF = ndims(data2D_mirror)-1;

% Each single loc (1-5)
for iLoc = 1:nLoc
    data2D_perComb(:, iLoc , :, :) = squeeze(data2D_mirror(:, iLoc, :, :));
    margOri_perComb{iLoc} = squeeze(mean(data2D_mirror(:, iLoc, :, :), margDimOri));
    margSF_perComb{iLoc} = squeeze(mean(data2D_mirror(:, iLoc, :, :), margDimSF));
    %     data2D_perLoc(:, iLoc, :, :) = data2D_comb{iLoc};% reorganize for saving processed 2D data
end

% HM (6)
data2D_perComb(:, 6, :, :) = squeeze(mean(data2D_mirror(:, [2,4], :, :), 2));
margOri_perComb{6} = squeeze(mean(mean(data2D_mirror(:, [2,4], :, :), 2), margDimOri));
margSF_perComb{6} = squeeze(mean(mean(data2D_mirror(:, [2,4], :, :), 2), margDimSF));

% VM (7)
data2D_perComb(:, 7, :, :) = squeeze(mean(data2D_mirror(:, [3,5], :, :), 2));
margOri_perComb{7} = squeeze(mean(mean(data2D_mirror(:, [3,5], :, :), 2), margDimOri));
margSF_perComb{7} = squeeze(mean(mean(data2D_mirror(:, [3,5], :, :), 2), margDimSF));

% Periphery (8)
data2D_perComb(:, 8, :, :) = squeeze(mean(data2D_mirror(:, 2:5, :, :), 2));
margOri_perComb{8} = squeeze(mean(mean(data2D_mirror(:, 2:5, :, :), 2), margDimOri));
margSF_perComb{8} = squeeze(mean(mean(data2D_mirror(:, 2:5, :, :), 2), margDimSF));

%% 8. fit perLoc and perComb
for ifeature = 1:2
    
    xaxis = xaxis_all{ifeature};
    nfilters = length(xaxis);
    nameModelParams = namesModelParams{ifeature};
    nparams = length(nameModelParams);
    
    % create empty containers
    marg_perComb_perF = nan(ntypes, ncomb, nfilters);
    margPred_perComb_perF = nan(ntypes, ncomb, nfilters); % PRS/ABS/BOTH, icomb, iline=1/2
    margParams_perComb_perF = nan(ntypes, ncomb, nparams);
    margR2_perComb_perF = nan(ntypes, ncomb);
    
    % model 
    imodelKernel = imodelKernel_perF(ifeature);
    ub = ub_all{imodelKernel};
    lb = lb_all{imodelKernel};
    
    switch ifeature
        case 1, marg_perComb_raw = margOri_perComb; 
        case 2, marg_perComb_raw = margSF_perComb; 
    end
    
    for itype = 1:ntypes
        for icomb = 1:8
            marg = marg_perComb_raw{icomb}(itype, :);
            switch ifeature
                case 1
                    [~, params_est, Rsquared, ~] = SX_sim08_fit(xaxis, marg, imodelKernel, ub, lb, fitMode);
                    margPred = predSFkernel(xaxis, imodelKernel, params_est, 0);
                case 2
                    [~, params_est, Rsquared, ~] = SX_sim08_fit(2.^xaxis, marg, imodelKernel, ub, lb, fitMode);
                    margPred = predSFkernel(2.^xaxis, imodelKernel, params_est, 0);
            end

            marg_perComb_perF(itype, icomb, :) = marg;
            margPred_perComb_perF(itype, icomb, :) = margPred;
            margParams_perComb_perF(itype, icomb, :) = params_est;
            margR2_perComb_perF(itype, icomb) = Rsquared;
        end
    end
    marg_perComb{ifeature} = marg_perComb_perF;
    margPred_perComb{ifeature} = margPred_perComb_perF;
    margParams_perComb{ifeature} = margParams_perComb_perF;
    margR2_perComb{ifeature} = margR2_perComb_perF;
end






