function [e3D, iPRS, resp, cst, respC, iPair] = fxn_getRespC(e3D_nonrand, indUnikPair_rand, iPair_nonrand, iPair_unik_nonrand, iPRS_nonrand, resp_nonrand, cst_nonrand)

% indRand = indRand_unikPair_train; % range: [1, ntrialsPerLoc/2], not repeated
iPair_unik = iPair_unik_nonrand(indUnikPair_rand);
nUNIK = length(iPair_unik);
indTrial = nan(nUNIK, 2);
for iu = 1:nUNIK
    indTrial(iu, :) = find(iPair_unik(iu) == iPair_nonrand);
end

e3D = e3D_nonrand(indTrial(:), :, :);
iPRS_2 = [iPRS_nonrand(indTrial(:, 1)), iPRS_nonrand(indTrial(:, 2))]; assert(sum(iPRS_2(:, 1) ~= iPRS_2(:, 2)) == 0)
resp_2 = [resp_nonrand(indTrial(:, 1)), resp_nonrand(indTrial(:, 2))];
iPRS = iPRS_nonrand(indTrial(:), :, :);
resp = resp_nonrand(indTrial(:), :, :);
cst = cst_nonrand(indTrial(:));
iPair = iPair_nonrand(indTrial(:));

%% get pA
respC = nan(nUNIK, 3); % pA_both/PRS/ABS (must be in this order)
for iu = 1:nUNIK
    respA = resp_2(iu, 1);
    respB = resp_2(iu, 2);
    assert(iPRS_2(iu, 1) == iPRS_2(iu, 2))
    if iPRS_2(iu, 1) == 1, respC(iu,:) = [respA==respB, respA==respB, nan]; % PRS
    else, respC(iu,:) = [respA==respB, nan,respA==respB]; % ABS
    end
end
