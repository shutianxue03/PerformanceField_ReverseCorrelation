


% INPUTS:
% energyFlag, energy_pp_perLoc/energy_noise_perLoc/energy_noiseMask_perLoc,
% energy_norm_perLoc/energy_norm_perLoc/energy_norm_perLoc
% ntypes, nfiltersSF, nfiltersOri,  nLoc, nbins_e, energy_perLoc, data_both

% INSIDE FXN
% checkEnerBins()
function [pYES_tgt, ebin_tgt] = SX_RC5_slope(dataMatrix, e3D, ntypes, nbins_e)
%% switch the type of energy used (patch, pure noise or masked noise)
% switch energyFlag
%     case 1 % energy is calculated by the presented patch:
%         energy_perLoc = energy2D_allT_perLoc; % 'allT' stands for 'all trials'
%         energy_norm_perLoc = energy2D_norm_allT_perLoc;
%     case 2 % energy is calculated by pure noise:
%         energy_perLoc = energy2D_noise_allT_perLoc;
%         energy_norm_perLoc = energy2D_noise_norm_allT_perLoc;
%     case 3 % energy is calculated by modified noise:
%         energy_perLoc = energy2D_noiseMask_allT_perLoc;
%         energy_norm_perLoc = energy2D_noiseMask_norm_allT_perLoc;
% end

%% ebin and pC at gabor SF and ORI
% use unnormalized energy; though the fxn shape won't change using normalized energy
[ntrials, nORI, nSF] = size(e3D); assert(size(dataMatrix, 1) == size(e3D, 1))
pYES_tgt = nan(ntypes, nbins_e);
ebin_tgt = pYES_tgt;

for itype = 1:ntypes
    switch itype
        case 1, iPRS = dataMatrix(:, 6) == 1;
        case 2, iPRS = dataMatrix(:, 6) == 0;
        case 3, iPRS = 1:ntrials;
    end
    [energy_allBins, pYES_allBins] = checkEnerBins(e3D(iPRS, ceil(nORI/2), ceil(nSF/2)), dataMatrix(iPRS, 9), nbins_e);
    pYES_tgt(itype, :) = pYES_allBins;
    ebin_tgt(itype, :) = energy_allBins;
    %         fprintf('Done iLoc = %d itype = %d\n', iLoc, itype)
end

% % plot
r = pYES_tgt;
e = ebin_tgt;
errType=1; ntrialsProp = 1;
% quickPlot5_slope

%% check the rHit/FA vs. energy for SELECTED ORI/SF filters
% checkSlopesAllFilters

%% estimate slope for each SF and ORI (use non-standardized energy)
% slopes2D_perComb = nan(2, nLoc8, nfiltersOri, nfiltersSF); % change to (nLoc, PRS/ABS, nfiltersOri, nfiltersSF)
% R2_2D_slope_perComb = slopes2D_perComb;
%
% for itype = 1:ntypes
%     switch itype
%         case 1, itrialStart = 1; itrialEnd = ntrials;
%         case 2, itrialStart = ntrials+1; itrialEnd = ntrials*2;
%         case 3, itrialStart = 1; itrialEnd = ntrials*2;
%     end
%
%     for ifilterOri = 1:nfiltersOri
%         for ifilterSF = 1:nfiltersSF
%             for iLoc = 1:nLoc8
%                 [pC_allBins, energy_allBins] = checkEnerBins(energy_perComb{iLoc}(itrialStart:itrialEnd, ifilterOri, ifilterSF), data_perComb{iLoc}(itrialStart:itrialEnd, 2), nbins_e);
%
%                 % get slopes
%                 [beta, ~] = polyfit(energy_allBins, pC_allBins, 1); % linear fitting
%                 % get slope
%                 slopes2D_perComb(itype, iLoc, ifilterOri, ifilterSF) = beta(1); % beta = [slope, intercept]
%                 % get R2
%                 yfit = polyval(beta, energy_allBins);
%                 R2_2D_slope_perComb(itype, iLoc, ifilterOri, ifilterSF) = getR2(pC_allBins, yfit);
%             end
%         end
%     end
% end


%% separability
% seprblity_slope_perComb = nan(2, nLoc8);
% for itype = 1:ntypes
%     for iLoc = 1:nLoc8
%         seprblity_slope_perComb(itype, iLoc) = getSeparability(squeeze(slopes2D_perComb(itype, iLoc, :,:)));
%     end
% end



