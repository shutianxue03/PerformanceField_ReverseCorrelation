function [e3D_rand, dataMtx_rand, respC] = fxn_getRespC_simp(indRand_unikPair, iPair_unik_nonrand, e3D_nonrand, dataMtx_nonrand)

% the simplied version of fxn_model/fxn_getRespC
% the input is just the dataMatrix

iPRS_nonrand = dataMtx_nonrand(:, 6);
iPair_nonrand = dataMtx_nonrand(:, 8);
resp_nonrand = dataMtx_nonrand(:, 9);

iPair_unik_rand = iPair_unik_nonrand(indRand_unikPair);
nUNIK = length(iPair_unik_rand);
indTrial = nan(2, nUNIK);
for iu = 1:nUNIK
    indTrial(:, iu) = find(iPair_unik_rand(iu) == iPair_nonrand);
end

e3D_rand = e3D_nonrand(indTrial(:), :, :);
dataMtx_rand = dataMtx_nonrand(indTrial(:), :);

iPRS_rand = [iPRS_nonrand(indTrial(1,:))'; iPRS_nonrand(indTrial(2,:))']; assert(sum(iPRS_rand(1, :) ~= iPRS_rand(2, :)) == 0)
resp_rand = [resp_nonrand(indTrial(1, :))'; resp_nonrand(indTrial(2, :))'];
% cst_rand = cst_nonrand(indTrial(:));

%% get pA
respC = nan(nUNIK, 3); % pA_both/PRS/ABS (must be in this order)
for iu = 1:nUNIK
    respA = resp_rand(1,iu);
    respB = resp_rand(2,iu);
    assert(iPRS_rand(1, iu) == iPRS_rand(2, iu))
    if iPRS_rand(1, iu) == 1, respC(iu,:) = [respA==respB, respA==respB, nan]; % PRS
    else, respC(iu,:) = [respA==respB, nan,respA==respB]; % ABS
    end
end
