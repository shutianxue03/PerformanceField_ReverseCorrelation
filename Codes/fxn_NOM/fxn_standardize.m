function [IV_PRS_norm, IV_ABS_norm] = fxn_standardize(IV_PRS, IV_ABS)
IV_mean = mean([IV_PRS, IV_ABS]);
IV_std = std([IV_PRS, IV_ABS]);
IV_PRS_norm = (IV_PRS - IV_mean)/IV_std;
IV_ABS_norm = (IV_ABS - IV_mean)/IV_std;
