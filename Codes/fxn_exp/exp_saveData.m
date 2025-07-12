

itrial = run.itrial;

record.iCuedLoc(iblock_current,itrial) = run.iCuedLoc;
record.answer(iblock_current,itrial) = run.answer;
record.correctness(iblock_current, itrial) = run.correctness;
record.tgtPrs(iblock_current, itrial) = run.tgtPrs;
record.RT(iblock_current,itrial) = run.RT;
record.phase_tgt(iblock_current,itrial) = run.phase_tgt;

if exp_mode ~= 1
    ib = iblock_current; % iblock_current is not record.iblock_current, which is iblock_current+1 already
    record.patch_allTrials = [];
    record.run_allTrials{ib, itrial} = run;
    record.iStair(iblock_current,itrial) = run.iStair;
    record.stairs_all{run.iCuedLoc,run.iStair} = run.stair;
else % record the repetition Ind to examine response consistency
    record.repInd(iblock_current, itrial) = run.repInd;
end