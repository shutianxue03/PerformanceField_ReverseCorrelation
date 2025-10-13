

function pred = fxn_getEstimates(params_est, IV_PRS, IV_ABS, plotFlag)
% function pred = fxn_getEstimates(params_est, IV_PRS, IV_ABS, ntrials_pA, normIV, plotFlag)
% do not care about the loc for now

% plot histogram
if plotFlag
    ModelPlot_hist(1, 1, IV_PRS, IV_ABS, params_est)
end

% predict metrics



