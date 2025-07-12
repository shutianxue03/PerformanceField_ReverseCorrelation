

% pC
pC_perType(1,ib) = nanmean(record.correctness(ib, boolean(record.tgtPrs(ib, :))));   % when gabor is PRS, i.e., rHit
pC_perType(2,ib) = nanmean(record.correctness(ib, boolean(1-record.tgtPrs(ib, :)))); % when gabor is ABS, i.e., rFA
pC_perType(3,ib) = nanmean(record.correctness(ib, :));    % all trials

assert(pC_perType(3,ib) - mean(pC_perType([1,2], ib))<eps) % assert that pC = (rHit+rCR)/2

% pA
for irep = repInd_all
    indRep = record.repInd(ib, :) == irep; % the pair of repeated trial
    ansRep = record.answer(ib, indRep); % the correctness of two responses
    respConsistency = ansRep(1) ==ansRep(2);
    
    iPrs = record.tgtPrs(ib, indRep); % whether this pair of trials is PRS/ABS
    if  iPrs(1) == 1
        pA_perType(1, ib-ib_start+1, irep == repInd_all) = respConsistency;
    elseif iPrs(1) == 0
        pA_perType(2, ib-ib_start+1, irep == repInd_all) = respConsistency;
    end
    pA_perType(3, ib-ib_start+1, irep == repInd_all) = respConsistency;
end

