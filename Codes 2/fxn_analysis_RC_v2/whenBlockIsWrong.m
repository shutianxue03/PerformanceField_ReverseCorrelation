

% step 1: load the last correct files, rename the record to record_c
% step 2: load the first wrong file, rename the record to record_w, adjust line 8-11
% of this script
% step 3: run this script 

% diff = -15;
subjName_c = 'RE';
ind_c = 119; % the ib supposed to be
ind_w = 120; % to locate the line in the wrong data matrix
loc_c = 2;
time_c = '20211105T1837'; 

%
fileName.formal = sprintf('Data/%s/%s_exp_B0%dL%d_%s', subjName_c, subjName_c, ind_c, loc_c, time_c);

record_c.iblock_current = ind_c + 1;
record_c.iCuedLoc(ind_c, :) = record_w.iCuedLoc(ind_w, :);
record_c.RT(ind_c, :) = record_w.RT(ind_w, :);
record_c.answer(ind_c, :) = record_w.answer(ind_w, :);
record_c.correctness(ind_c, :) = record_w.correctness(ind_w, :);
record_c.tgtPrs(ind_c, :) = record_w.tgtPrs(ind_w, :);
record_c.phase_tgt(ind_c, :) = record_w.phase_tgt(ind_w, :);
record_c.repInd(ind_c, :) = record_w.repInd(ind_w, :);

record_c.thresh_exp = record_w.thresh_exp;
record_c.accPerBlock = record_w.accPerBlock;
record_c.dprime = record_w.dprime;
record_c.criterion = record_w.criterion;
record_c.run_allTrials = record_w.run_allTrials;
record_c.duration = record_w.duration;
record_c.fileName = fileName;

record = record_c; 

save(fileName.formal, 'exp_mode', 'fileName', 'params', 'record')



