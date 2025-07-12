function run = exp_checkResp(params, run, giveFeedback)

stim = params.stim;
time = params.time;
soundParams = params.soundParams;

if run.key == stim.keysCodes(1), run.answer = 1; else, run.answer = 0; end
% run.answer = 1; % mute this line if you want to make response manually
if run.tgtPrs == run.answer
    run.correctness = 1;
    if giveFeedback
        makeBeep(soundParams, soundParams.freqCorrect, time.feedback)
    else,WaitSecs(time.feedback);
    end
else
    run.correctness = 0;
    if giveFeedback
        makeBeep(soundParams, soundParams.freqWrong, time.feedback)
    else,WaitSecs(time.feedback);
    end
end
