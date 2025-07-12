function run = SX9_initRun(itrial, iblock, design, record, exp_mode)
% Note that the matrix 'design.ntrial_ind' contains parameters of all
% trials of all blocks.

run.iblock = iblock;
run.itrial = itrial;
realTrial = (iblock-1)*design.nTrialsPerBlock + itrial;

if exp_mode == 1 % also extract the created stimuli
    run_extracted = record.run_allTrials{run.itrial};
    run.tgtPrs = run_extracted.tgtPrs;         % whether a gabor was presented at the tgt loc
%     run.tgtPrs_all = run_extracted.tgtPrs_all;% whether a gabor was presented at all loc
    run.iCuedLoc = run_extracted.iCuedLoc;
    run.stimPatch_all = run_extracted.stimPatch_all;    % all the stimuli presented
    run.noisePatch_all = run_extracted.noisePatch_all; % all the noise patches used to make a stimulus
    run.phase_tgt = run_extracted.phase_tgt;
    run.repInd = run_extracted.repInd;
    run.iStair = nan;
    run.stair = nan;
else
    run.iCuedLoc = design.trialInd(realTrial,1);   % location of the cue
    run.iStair = design.trialInd(realTrial,2);         % which staircase at this location (1:2)
    run.tgtPrs = design.trialInd(realTrial,3);        % whether the cued stim is a target or not (1,0)
%     run.tgtPrs_all = design.
    run.stair = record.stairs_all{run.iCuedLoc, run.iStair};
    run = exp_drawStim(exp_mode, record.params, run, record);
end


run.saveBlock = 0;
run.thresh_exp = record.thresh_exp;


