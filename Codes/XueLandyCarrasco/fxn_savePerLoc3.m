% function fxn_savePerLoc

% metrics
similarityAfterMirroring_3allSubj(isubj, iB_start:iB_end, iiLoc3, :) = similarityAfterMirroring_allB(:, :, iiLoc2);
sep_3allSubj(isubj, iB_start:iB_end, iiLoc3, :) = sep_allB(:, iiLoc2, :);
% marg
margORI_3allSubj(isubj, iB_start:iB_end, iiLoc3, :, :) = margORI_allB(:, iiLoc2, :, :);
margSF_3allSubj(isubj, iB_start:iB_end, iiLoc3, :, :) = margSF_allB(:, iiLoc2, :, :);
% pred
margPredORI_3allSubj(isubj, iB_start:iB_end, iiLoc3, :, :) = margPredORI_allB(:, iiLoc2, :, :);
margPredSF_3allSubj(isubj, iB_start:iB_end, iiLoc3, :, :) = margPredSF_allB(:, iiLoc2, :, :);
margR2ORI_3allSubj(isubj, iB_start:iB_end, iiLoc3, :, :) = margR2ORI_allB(:, iiLoc2, :);
margR2SF_3allSubj(isubj, iB_start:iB_end, iiLoc3, :, :) = margR2SF_allB(:, iiLoc2, :);
% estimated params
margParamsORI_3allSubj(isubj, iB_start:iB_end, iiLoc3, :, :) = margParamsORI_allB(:, iiLoc2, :, :);
margParamsSF_3allSubj(isubj, iB_start:iB_end, iiLoc3, :, :) = margParamsSF_allB(:, iiLoc2, :, :);
% tuning characteristics
margTuningC_ORI_3allSubj(isubj, iB_start:iB_end, iiLoc3, :, :) = margTuningC_ORI_allB(:, iiLoc2, :, :);
margTuningC_SF_3allSubj(isubj, iB_start:iB_end, iiLoc3, :, :) = margTuningC_SF_allB(:, iiLoc2, :, :);

if flag_ABSprefORI
    margTuningC_ORI_3allSubj(:, :, :, :, 1) = abs(margTuningC_ORI_3allSubj(:, :, :, :, 1));
end