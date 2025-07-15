% function fxn_savePerLoc
similarityAfterMirroring_4allSubj(isubj, :, iiLoc4, :) = similarityAfterMirroring_allB(:, :, iiLoc2);
% kernels2D_4allSubj(isubj, :, iiLoc4, :, :, :) = kernels2D_allB(:, iiLoc2, :, :, :);
sep_4allSubj(isubj, :, iiLoc4, :) = sep_allB(:, iiLoc2, :);
% marg
margORI_4allSubj(isubj, :, iiLoc4, :, :) = margORI_allB(:, iiLoc2, :, :);
margSF_4allSubj(isubj, :, iiLoc4, :, :) = margSF_allB(:, iiLoc2, :, :);
% pred
margPredORI_4allSubj(isubj, :, iiLoc4, :, :) = margPredORI_allB(:, iiLoc2, :, :);
margPredSF_4allSubj(isubj, :, iiLoc4, :, :) = margPredSF_allB(:, iiLoc2, :, :);
% % estimated params
margParamsORI_4allSubj(isubj, :, iiLoc4, :, :) = margParamsORI_allB(:, iiLoc2, :, :);
margParamsSF_4allSubj(isubj, :, iiLoc4, :, :) = margParamsSF_allB(:, iiLoc2, :, :);
% tuning characteristics
margTuningC_ORI_4allSubj(isubj, :, iiLoc4, :, :) = margTuningC_ORI_allB(:, iiLoc2, :, :);
margTuningC_SF_4allSubj(isubj, :, iiLoc4, :, :) = margTuningC_SF_allB(:, iiLoc2, :, :);

