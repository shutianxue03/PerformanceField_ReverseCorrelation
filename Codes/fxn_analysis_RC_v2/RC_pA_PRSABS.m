% 
% 
% 
% % repInd_all = unique(record.repInd(1,:));
% repInd_all = 1:50;
% ntrialsRep = length(repInd_all);
% 
% % just calculate a rough pA without considering signal presence
% pA_ = nan(ntypes, nblocks, ntrialsRep);
% pC_ = nan(ntypes, nblocks);
% 
% for ib = ib_start:ib_end
% 
% %     % pC
% %     pC_(1,ib) = nanmean(record.correctness(ib, boolean(record.tgtPrs(ib, :))));   % gabor-PRS
% %     pC_(2,ib) = nanmean(record.correctness(ib, boolean(1-record.tgtPrs(ib, :)))); % gabor-ABS   
% %     pC_(3,ib) = nanmean(record.correctness(ib, :));    % both
% %     assert(pC_(3,ib) - mean(pC_([1,2], ib))<eps) % it is supposed to be pC = (rHit+rCR)/2
% %     
% %     % pA
% %     for irep = repInd_all
% %         indRep = record.repInd(ib, :) == irep; % the pair of repeated trial
% %         ansRep = record.answer(ib, indRep); % the correctness of two responses
% %         respConsistency = ansRep(1) ==ansRep(2);
% %         
% %         iPrs = record.tgtPrs(ib, indRep); % whether this pair of trials is PRS/ABS
% %         if  iPrs(1) == 1
% %             pA_(1, ib-ib_start+1, irep == repInd_all) = respConsistency;
% %         elseif iPrs(1) == 0
% %             pA_(2, ib-ib_start+1, irep == repInd_all) = respConsistency;
% %         end
% %         
% %         pA_(3, ib-ib_start+1, irep == repInd_all) = respConsistency;
% %         
% %     end
% end
% 
% % pA = nanmean(pA_, 3);
% % 
% % pA_perSess_perLoc = nan(ntypes, nblocks/nLoc8, nLoc8);
% % pC_perSess_perLoc = pA_perSess_perLoc;
% % for iLoc = 1:nLoc8
% %     pA_perSess_perLoc(:, :, iLoc) = pA(:, iCuedLoc==iLoc); 
% %     pC_perSess_perLoc(:, :, iLoc) = pC_(:, iCuedLoc==iLoc); 
% % end
% 
