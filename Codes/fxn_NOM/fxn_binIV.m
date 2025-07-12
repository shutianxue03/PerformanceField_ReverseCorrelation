function [nTrials_allBins, ind] = fxn_binIV(IV, nBins)
[nTrials_allBins, ~, ind] = histcounts(IV, nBins);
end