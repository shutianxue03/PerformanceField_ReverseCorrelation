

%% calculate pA, regardless of signal existence
repInd_all = unique(record.repInd(1,:));
ntrialsRep = length(repInd_all);

% just calculate a rough pA without considering signal presence
pA_ = nan(nblocks, ntrialsRep);

for ib = ib_start:ib_end
    for repInd = repInd_all
        ansRep = record.answer(ib, record.repInd(ib, :) == repInd); % the correctness of two responses
        pA_(ib-ib_start+1, repInd == repInd_all) = (ansRep(1) ==ansRep(2));
    end
end

pA = mean(pA_, 2);

pA_perSess_perLoc = nan(nblocks/nLoc, nLoc);
for iLoc = 1:nLoc, pA_perSess_perLoc(:, iLoc) = pA(iCuedLoc==iLoc); end
pA_perSess_perLoc_ave = mean(pA_perSess_perLoc);
pA_perSess_perLoc_sd = std(pA_perSess_perLoc)/sqrt(nblocks/nLoc);


%% calculate pA, considering signal existence
nresp2_abs = nan(nblocks,3);
nresp2_prs = nresp2_abs;

for ib = ib_start:ib_end
    ansPrs = nan(ntrialsRep/2, 2);
    ansAbs = ansPrs;
    ansPrsInd = 1;
    ansAbsInd = 1;
    
    for repInd = repInd_all
        ansPrs_ = record.answer(ib, record.repInd(ib, :) == repInd & record.tgtPrs(ib, :) == 1);
        ansAbs_ = record.answer(ib, record.repInd(ib, :) == repInd & record.tgtPrs(ib, :) == 0);
        
        if isempty(ansPrs_) % the current pair are signal-absent
            ansAbs(ansAbsInd, :) = ansAbs_;
            ansAbsInd = ansAbsInd+1;
        else
            ansPrs(ansPrsInd, :) = ansPrs_;
            ansPrsInd = ansPrsInd+1;
        end
    end

    nresp2_prs(ib-ib_start+1, 1) = sum(ansPrs(:,1) == 0 & ansPrs(:,2) == 0); % both respond prs given signal-prs
    nresp2_prs(ib-ib_start+1, 2) = sum(ansPrs(:,1) ~= ansPrs(:,2)); % inconsistent response
    nresp2_prs(ib-ib_start+1, 3) = sum(ansPrs(:,1) == 1 & ansPrs(:,2) == 1); % both respond abs given signal-prs
    
    nresp2_abs(ib-ib_start+1, 1) = sum(ansAbs(:,1) == 0 & ansAbs(:,2) == 0); % both respond PRS given signal-abs
    nresp2_abs(ib-ib_start+1, 2) = sum(ansAbs(:,1) ~= ansAbs(:,2)); % inconsistent response
    nresp2_abs(ib-ib_start+1, 3) = sum(ansAbs(:,1) == 1 & ansAbs(:,2) == 1); % both respond abs given signal-abs
    
end

pA_prs = sum(nresp2_prs)./(ntrialsRep/2*nblocks);
pA_abs = sum(nresp2_abs)./(ntrialsRep/2*nblocks);

% confirm that the resp consistency is the same
pA2 = (nresp2_prs(:,1) + nresp2_prs(:,3) + nresp2_abs(:,1) + nresp2_abs(:,3))./sum(nresp2_prs+nresp2_abs,2);
assert(pA2(1) == pA(1) & pA2(2) == pA(2))

nresp2_prs_byLoc = nan(3, nLoc);
nresp2_abs_byLoc = nresp2_prs_byLoc;
for iLoc = 1:nLoc
    nresp2_prs_byLoc(:, iLoc) = sum(nresp2_prs(iCuedLoc==iLoc, :));
    nresp2_abs_byLoc(:, iLoc) = sum(nresp2_abs(iCuedLoc==iLoc, :));
end

pA_prs_perLoc =  nresp2_prs_byLoc./(ntrialsRep/2*nblocks/nLoc);
pA_abs_perLoc =  nresp2_abs_byLoc./(ntrialsRep/2*nblocks/nLoc);








