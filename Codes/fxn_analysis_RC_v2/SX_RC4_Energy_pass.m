
% INPUTS:
% nLoc, patch_both, noise_both, noise_mask_both, cst_both, nAllTrials,
% nfiltersOri, nfiltersSF, filter_sin, filter_cos, filtersSF_all, filtersOri_all

% INSIDE FXN
% SX_sim04_computeEnergy()
% quickPlot_energy()
function [energy2D_passA, energy2D_passB] = SX_RC4_Energy_pass(mask, patch_both, data_both, filter_sin, filter_cos)

[nfiltersOri, nfiltersSF] = size(filter_sin);
plot4parts = 0;     % 1 = plot the stim and two filters
nLoc8 = 8;

e = cell(nLoc8, 1);
e_norm = e;
phase = e;

%% collect energy by combining location
iLoc = 1;
% for iLoc = 1:nLoc8
% extract data
patch_perLoc = patch_both{iLoc}; % name it 'patch_' since patch is a function
cst_perLoc = data_both{iLoc, 5}; % 5th column is CST
passInd_perLoc = data_both{iLoc, 6}; % 6th: pass ind (1 = pass A, 2 = pass B)
lumiBG = patch_perLoc{1}(1,1);
nAllTrials = length(patch_perLoc);

% empty containers

energy2D_allT = nan(nAllTrials, nfiltersOri, nfiltersSF);
% stimPhase2D_allT = energy2D_allT;

% loop to calculate energy
for ifilterOri = 1:nfiltersOri
    for ifilterSF = 1:nfiltersSF
        filter_sin_ = filter_sin{ifilterOri, ifilterSF};
        filter_cos_ = filter_cos{ifilterOri, ifilterSF};
        parfor itrial = 1:nAllTrials
            [a,b] = SX_sim04_computeEnergy(mask, patch_perLoc{itrial} , lumiBG, filter_sin_, filter_cos_, plot4parts);
            energy2D_allT(itrial, ifilterOri, ifilterSF) = a;
        end
    end
    fprintf('Loc #%d/%d: ORI #%d/%d SF #%d/%d DONE.\n', iLoc, nLoc8, ifilterOri, nfiltersOri, ifilterSF, nfiltersSF)
end

% save
energy2D_passA = energy2D_allT(passInd_perLoc == 1, :, :);
energy2D_passB = energy2D_allT(passInd_perLoc == 2, :, :);




