function run = test_constructRun(ntrial,testStruct)
run.ntrial = ntrial;
run.nBlock = ceil(ntrial/testStruct.nTrialPerBlock);
run.nTargetLoc = testStruct.ntrial_ind(ntrial,1); % location of the target
run.nCueLoc = testStruct.ntrial_ind(ntrial,2);  % location of the cue
run.noiseSF = testStruct.ntrial_ind(ntrial,3); % SF of the noise
run.noiseOri = testStruct.ntrial_ind(ntrial,4); % orientation of the noise
run.trueAnswer = (run.nTargetLoc == run.nCueLoc);  % 0 = incorrect, 1 = correct
run.saveBlock = 0;