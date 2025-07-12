

% run_allTrials: col cell [ntrialsPerBlock x 1]
%     - each cell contains a struct 'run', each contains stimPatch,
%     noisePatch, iCuedLoc, contrast, iTratget, phase_tgt
design.nrep = 2;
ntrialsPerRep = design.nTrialsPerBlock/design.nrep;
run = [];

run_all1_nonrand = cell(ntrialsPerRep, 1); 
run.iCuedLoc = design.loc_list(idx_block);

for n = [1,0] % present and absent
    run.tgtPrs = n;
    if mod(ntrialsPerRep,2) ~= 0, error('ntrials is not multiple of 4!! '), end
    for ii = 1:ntrialsPerRep/2
        ind = ii+(1-n)*ntrialsPerRep/2;
        run_all1_nonrand{ind} = exp_drawStim(exp_mode, params, run, record);
        run_all1_nonrand{ind}.repInd = ind;
    end
end

% randomize sequence and combine two reps
run_all1_ind = randperm(ntrialsPerRep);
run_all1 = run_all1_nonrand(run_all1_ind); % randomize the sequence of the 1st rep
run_all2_ind = run_all1_ind(randperm(ntrialsPerRep));
run_all2 = run_all1(run_all2_ind, :); % copy the first rep and randomize the sequence
run_allTrials = [run_all1; run_all2];


