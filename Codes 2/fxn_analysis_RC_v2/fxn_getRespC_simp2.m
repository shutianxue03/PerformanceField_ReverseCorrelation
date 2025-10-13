function respC = fxn_getRespC_simp2(iPair, iPRS, resp)

% a simplified version to get resp consistency matrix only
iPair_unik = unique(iPair);
nUNIK = length(iPair_unik);
indTrial = nan(nUNIK, 2);
for iu = 1:nUNIK
    indTrial(iu, :) = find(iPair_unik(iu) == iPair);
end

iPRS_2 = [iPRS(indTrial(:, 1)), iPRS(indTrial(:, 2))]; assert(sum(iPRS_2(:, 1) ~= iPRS_2(:, 2)) == 0)
resp_2 = [resp(indTrial(:, 1)), resp(indTrial(:, 2))];
iPRS = iPRS(indTrial(:), :, :);
resp = resp(indTrial(:), :, :);
iPair = iPair(indTrial(:));

%% get resp C
respC = nan(nUNIK, 3); % pA_both/PRS/ABS (must be in this order)
for iu = 1:nUNIK
    respA = resp_2(iu, 1);
    respB = resp_2(iu, 2);
    assert(iPRS_2(iu, 1) == iPRS_2(iu, 2))
    if iPRS_2(iu, 1) == 1, respC(iu,:) = [respA==respB, respA==respB, nan]; % PRS
    else, respC(iu,:) = [respA==respB, nan,respA==respB]; % ABS
    end
end
